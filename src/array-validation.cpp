#include <Rcpp.h>
#include <unordered_map>
#include <unordered_set>
using namespace Rcpp;

namespace {
struct ArrayValidator {
  Function mono, sub, link, anomer, floating;
  std::unordered_map<SEXP, bool> monos, subs, links, anomers, floats;
  ArrayValidator(Function m, Function s, Function l, Function a, Function f)
      : mono(m), sub(s), link(l), anomer(a), floating(f) {}
  bool cached(SEXP s, Function& fn, std::unordered_map<SEXP, bool>& cache) {
    // Input CHARSXPs stay rooted in the records for this batch. Distinct
    // encodings may miss the cache but still use the same R validator.
    SEXP key = s;
    auto found = cache.find(key);
    if (found != cache.end()) return found->second;
    bool ok = as<bool>(fn(CharacterVector::create(s)));
    cache.emplace(key, ok);
    return ok;
  }
  bool chars(SEXP x, R_xlen_t n = -1) {
    if (TYPEOF(x) != STRSXP || (n >= 0 && XLENGTH(x) != n)) return false;
    for (R_xlen_t i = 0; i < XLENGTH(x); ++i)
      if (STRING_ELT(x, i) == NA_STRING) return false;
    return true;
  }
  bool indices(SEXP x, int n, bool unique = false) {
    if (TYPEOF(x) != INTSXP) return false;
    std::unordered_set<int> seen;
    for (R_xlen_t i = 0; i < XLENGTH(x); ++i) {
      int id = INTEGER(x)[i];
      if (id < 1 || id > n || (unique && !seen.insert(id).second)) return false;
    }
    return true;
  }
  bool fields(SEXP x, std::initializer_list<const char*> allowed,
              int required, std::vector<SEXP>& out) {
    if (TYPEOF(x) != VECSXP) return false;
    SEXP names = Rf_getAttrib(x, R_NamesSymbol);
    if (TYPEOF(names) != STRSXP || XLENGTH(names) != XLENGTH(x)) return false;
    out.assign(allowed.size(), R_NilValue);
    std::vector<bool> seen(allowed.size(), false);
    for (R_xlen_t i = 0; i < XLENGTH(x); ++i) {
      if (STRING_ELT(names, i) == NA_STRING) return false;
      const char* name = CHAR(STRING_ELT(names, i));
      auto at = std::find_if(allowed.begin(), allowed.end(),
        [&](const char* value) { return std::strcmp(value, name) == 0; });
      if (at == allowed.end()) return false;
      size_t j = at - allowed.begin();
      if (seen[j]) return false;
      seen[j] = true; out[j] = VECTOR_ELT(x, i);
    }
    for (int j = 0; j < required; ++j) if (!seen[j]) return false;
    return true;
  }
  int check(SEXP record) {
    if (Rf_isNull(record)) return 1;
    std::vector<SEXP> a;
    if (!fields(record, {"mono", "sub", "edges", "linkage", "anomer", "alditol",
                         "floating_parts", "floating_substituents"}, 6, a)) return 0;
    if (!chars(a[0]) || XLENGTH(a[0]) < 1 || XLENGTH(a[0]) > INT_MAX) return 0;
    int n = XLENGTH(a[0]);
    if (!chars(a[1], n) || !indices(a[2], n) || XLENGTH(a[2]) % 2 ||
        !chars(a[3], XLENGTH(a[2])/2) || !chars(a[4], 1) ||
        TYPEOF(a[5]) != LGLSXP || XLENGTH(a[5]) != 1 || LOGICAL(a[5])[0] == NA_LOGICAL)
      return 0;
    std::vector<SEXP> part_links, float_subs;
    for (int k = 6; k <= 7; ++k) {
      if (Rf_isNull(a[k])) continue;
      if (TYPEOF(a[k]) != VECSXP) return 0;
      for (R_xlen_t i = 0; i < XLENGTH(a[k]); ++i) {
        std::vector<SEXP> p;
        if (k == 6) {
          if (!fields(VECTOR_ELT(a[k], i), {"root", "nodes", "linkage", "parents"}, 4, p) ||
              !indices(p[0], n) || XLENGTH(p[0]) != 1 || !indices(p[1], n, true) ||
              !chars(p[2], 1) || !indices(p[3], n, true)) return 0;
          part_links.push_back(STRING_ELT(p[2], 0));
        } else {
          if (!fields(VECTOR_ELT(a[k], i), {"substituent", "parents"}, 2, p) ||
              !chars(p[0], 1) || !indices(p[1], n, true)) return 0;
          float_subs.push_back(STRING_ELT(p[0], 0));
        }
      }
    }
    for (int i = 0; i < n; ++i) if (!cached(STRING_ELT(a[0], i), mono, monos)) return 2;
    for (int i = 0; i < n; ++i) if (!cached(STRING_ELT(a[1], i), sub, subs)) return 3;
    for (R_xlen_t i = 0; i < XLENGTH(a[3]); ++i)
      if (!cached(STRING_ELT(a[3], i), link, links)) return 4;
    if (!cached(STRING_ELT(a[4], 0), anomer, anomers)) return 5;
    for (SEXP p : part_links) if (!cached(p, link, links)) return 6;
    for (SEXP p : float_subs) if (!cached(p, floating, floats)) return 7;
    return 1;
  }
};
}

// Validate shared chemical tokens once per batch using the established R rules.
// Unusual schemas retain the reference validator and its diagnostics.
// [[Rcpp::export(name = ".compact_check_arrays_native")]]
IntegerVector compact_check_arrays_native(List records, Function mono, Function sub,
                                         Function link, Function anomer, Function floating) {
  ArrayValidator validator(mono, sub, link, anomer, floating);
  IntegerVector result(records.size());
  for (R_xlen_t i = 0; i < records.size(); ++i) {
    if (i % 128 == 0) checkUserInterrupt();
    result[i] = validator.check(records[i]);
  }
  return result;
}
