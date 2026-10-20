#include "compact-progress.h"
#include <cstdint>
#include <unordered_map>
using namespace Rcpp;

namespace {
// igraph's version-4 serialized representation stores zero-based endpoints and
// attributes in R vectors. Its C adjacency index is reconstructed by igraph on
// demand, just as after readRDS(). Never copy a live C pointer between graphs.
// A template obtained from the installed igraph must match the entire schema;
// otherwise callers fall back to its public constructors and attribute setters.
bool supported_template(SEXP x) {
  if (TYPEOF(x) != VECSXP || Rf_xlength(x) != 10 || !Rf_inherits(x, "igraph")) return false;
  List g(x);
  if (TYPEOF(g[0]) != REALSXP || Rf_xlength(g[0]) != 1 || as<double>(g[0]) != 0 ||
      TYPEOF(g[1]) != LGLSXP || Rf_xlength(g[1]) != 1 || as<int>(g[1]) != TRUE) return false;
  for (int i = 2; i < 4; ++i)
    if (TYPEOF(g[i]) != REALSXP || Rf_xlength(g[i]) != 0) return false;
  for (int i = 4; i < 8; ++i) if (!Rf_isNull(g[i])) return false;
  if (TYPEOF(g[8]) != VECSXP || Rf_xlength(g[8]) != 4 || TYPEOF(g[9]) != ENVSXP) return false;
  List attrs = g[8];
  if (TYPEOF(attrs[0]) != REALSXP || Rf_xlength(attrs[0]) != 3) return false;
  NumericVector flags = attrs[0];
  if (flags[0] != 1 || flags[1] != 0 || flags[2] != 1) return false;
  for (int i = 1; i < 4; ++i)
    if (TYPEOF(attrs[i]) != VECSXP || Rf_xlength(attrs[i]) != 0) return false;
  Environment env = g[9];
  if (!env.exists(".__igraph_version__.") || !env.exists("myid") || !env.exists("igraph")) return false;
  SEXP version = env[".__igraph_version__."];
  SEXP id = env["myid"];
  return TYPEOF(version) == INTSXP && Rf_xlength(version) == 1 && INTEGER(version)[0] == 4 &&
    TYPEOF(id) == STRSXP && Rf_xlength(id) == 1 && STRING_ELT(id, 0) != NA_STRING &&
    LENGTH(STRING_ELT(id, 0)) > 0;
}

List make_graph(const IntegerVector& edges, int n, bool directed, List attrs, List model) {
  static uint64_t serial = 0;
  if (edges.size() % 2 || n < 0) stop("Invalid compact graph dimensions.");
  NumericVector from(edges.size() / 2), to(edges.size() / 2);
  for (int i = 0; i < from.size(); ++i) {
    if (edges[2*i] < 1 || edges[2*i] > n || edges[2*i+1] < 1 || edges[2*i+1] > n)
      stop("Invalid compact graph endpoint.");
    from[i] = edges[2*i] - 1;
    to[i] = edges[2*i+1] - 1;
  }
  Environment source = model[9];
  Environment env(R_NewEnv(R_EmptyEnv, FALSE, 0));
  // The prefix is an igraph-generated UUID, and the counter is process-local.
  // Graph IDs are opaque strings; each graph has its own ID and environment.
  env["myid"] = as<std::string>(source["myid"]) + ":" + std::to_string(++serial);
  env[".__igraph_version__."] = static_cast<SEXP>(source[".__igraph_version__."]);
  env["igraph"] = R_NilValue;
  List graph = List::create((double)n, directed, from, to, R_NilValue,
    R_NilValue, R_NilValue, R_NilValue, attrs, env);
  graph.attr("class") = model.attr("class");
  return graph;
}
}

// [[Rcpp::export(name = ".compact_build_graphs_native")]]
SEXP compact_build_graphs_native(List records, SEXP model, SEXP progress = R_NilValue) {
  if (!supported_template(model)) return R_NilValue;
  List template_graph(model), template_attrs = template_graph[8];
  CompactProgress reporter(progress);
  List out(records.size());
  std::unordered_map<int, CharacterVector> node_names;
  for (int i = 0; i < records.size(); ++i) {
    if (i % 128 == 0) { checkUserInterrupt(); reporter.update(i); }
    List a = records[i];
    CharacterVector mono = a["mono"], subs = a["sub"], linkage = a["linkage"];
    IntegerVector edges = a["edges"];
    if (subs.size() != mono.size() || linkage.size() * 2 != edges.size())
      stop("Invalid compact graph attribute lengths.");
    auto entry = node_names.find(mono.size());
    if (entry == node_names.end()) {
      CharacterVector names(mono.size());
      for (int j = 0; j < mono.size(); ++j) names[j] = std::to_string(j+1);
      entry = node_names.emplace(mono.size(), names).first;
    }
    CharacterVector names = entry->second;
    List graph_attrs = List::create(_["anomer"] = a["anomer"], _["alditol"] = a["alditol"]);
    if (a.containsElementNamed("floating_parts")) graph_attrs["floating_parts"] = a["floating_parts"];
    if (a.containsElementNamed("floating_substituents")) graph_attrs["floating_substituents"] = a["floating_substituents"];
    NumericVector flags = clone(as<NumericVector>(template_attrs[0]));
    List attrs = List::create(flags, graph_attrs,
      List::create(_["name"] = names, _["mono"] = mono, _["sub"] = subs),
      List::create(_["linkage"] = linkage));
    out[i] = make_graph(edges, mono.size(), true, attrs, template_graph);
  }
  reporter.update(records.size(), true);
  return out;
}

// [[Rcpp::export(name = ".compact_blank_graph_native")]]
SEXP compact_blank_graph_native(IntegerVector edges, int n, bool directed, SEXP model) {
  if (!supported_template(model)) return R_NilValue;
  List template_graph(model);
  List attrs = clone(as<List>(template_graph[8]));
  return make_graph(edges, n, directed, attrs, template_graph);
}
