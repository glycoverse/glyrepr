#include <Rcpp.h>
#include <unordered_map>
#include <unordered_set>
using namespace Rcpp;

// [[Rcpp::export(name = ".graph_counts_native")]]
SEXP graph_counts_native(List records, CharacterVector components,
                         CharacterVector targets, CharacterVector generic,
                         bool count, bool uncertain) {
  std::unordered_map<std::string, int> index;
  std::unordered_set<std::string> selected, ambiguous;
  for (int j = 0; j < components.size(); ++j)
    index.emplace(as<std::string>(components[j]), j);
  for (auto value : targets) selected.insert(as<std::string>(value));
  for (auto value : generic) ambiguous.insert(as<std::string>(value));
  List compositions(records.size());
  IntegerVector totals(records.size());
  for (R_xlen_t i = 0; i < records.size(); ++i) {
    if (i % 128 == 0) checkUserInterrupt();
    List record = records[i];
    CharacterVector residues = record[0], subs = record[1];
    std::vector<int> counts(components.size(), 0);
    bool unknown = false;
    int total = 0;
    auto add = [&](const std::string& value) {
      auto found = index.find(value);
      if (found == index.end()) return false;
      ++counts[found->second];
      if (selected.count(value)) ++total;
      return true;
    };
    for (auto value : residues) {
      std::string residue = as<std::string>(value);
      if (!add(residue)) return R_NilValue;
      unknown = unknown || ambiguous.count(residue);
    }
    for (auto value : subs) {
      if (CharacterVector::is_na(value)) continue;
      std::string text = as<std::string>(value);
      size_t start = 0;
      while (start < text.size()) {
        size_t end = text.find(',', start);
        if (end == std::string::npos) end = text.size();
        if (end > start) {
          size_t letters = end;
          while (letters > start && ((text[letters - 1] >= 'A' && text[letters - 1] <= 'Z') ||
                 (text[letters - 1] >= 'a' && text[letters - 1] <= 'z'))) --letters;
          if (letters == end || !add(text.substr(letters, end - letters))) return R_NilValue;
        }
        start = end + 1;
      }
    }
    if (count) {
      totals[i] = residues.size() == 0 || (uncertain && unknown) ? NA_INTEGER : total;
    } else {
      IntegerVector values;
      CharacterVector names;
      for (int j = 0; j < components.size(); ++j) {
        if (counts[j]) { values.push_back(counts[j]); names.push_back(components[j]); }
      }
      values.attr("names") = names;
      compositions[i] = values;
    }
  }
  if (count) return totals;
  return compositions;
}
