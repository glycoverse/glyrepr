#include "compact-progress.h"
#include "compact-core.h"
#include <cstdint>
#include <Rversion.h>
#include <unordered_map>
using namespace Rcpp;

namespace glyrepr_compact {
Forest import_graph_data(SEXP mono, SEXP sub, SEXP anomer, SEXP alditol,
                         SEXP edge_data, SEXP link_data, SEXP parts, SEXP subs,
                         bool preserve_parent_order, SEXP to);
void transform_graph(Forest& f, const std::string& operation,
                     const std::unordered_map<std::string, std::string>& mapping);
}

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

List graph_from_endpoints(NumericVector from, NumericVector to, int n, bool directed,
                         List attrs, List model) {
  static uint64_t serial = 0;
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

List make_graph(const IntegerVector& edges, int n, bool directed, List attrs, List model) {
  if (edges.size() % 2 || n < 0) stop("Invalid compact graph dimensions.");
  NumericVector from(edges.size() / 2), to(edges.size() / 2);
  for (int i = 0; i < from.size(); ++i) {
    if (edges[2*i] < 1 || edges[2*i] > n || edges[2*i+1] < 1 || edges[2*i+1] > n)
      stop("Invalid compact graph endpoint.");
    from[i] = edges[2*i] - 1;
    to[i] = edges[2*i+1] - 1;
  }
  return graph_from_endpoints(from, to, n, directed, attrs, model);
}

bool no_attributes(SEXP x) {
#if R_VERSION >= R_Version(4, 6, 0)
  return !ANY_ATTRIB(x);
#else
  return ATTRIB(x) == R_NilValue;
#endif
}

bool plain(SEXP x, SEXPTYPE type, R_xlen_t size) {
  return TYPEOF(x) == type && Rf_xlength(x) == size && no_attributes(x);
}

bool strings(SEXP x, R_xlen_t size) {
  if (!plain(x, STRSXP, size)) return false;
  for (R_xlen_t i = 0; i < size; ++i)
    if (STRING_ELT(x, i) == NA_STRING) return false;
  return true;
}

bool fields(SEXP x, std::initializer_list<const char*> expected) {
  if (TYPEOF(x) != VECSXP || Rf_xlength(x) != (R_xlen_t)expected.size()) return false;
  SEXP names = Rf_getAttrib(x, R_NamesSymbol);
  if (TYPEOF(names) != STRSXP || Rf_xlength(names) != Rf_xlength(x)) return false;
  int i = 0;
  for (const char* name : expected)
    if (STRING_ELT(names, i) == NA_STRING ||
        std::strcmp(CHAR(STRING_ELT(names, i++)), name)) return false;
  return true;
}

bool indices(SEXP x, int n, bool nonempty = false) {
  if (TYPEOF(x) != INTSXP || !no_attributes(x) ||
      (nonempty && Rf_xlength(x) == 0)) return false;
  for (int id : IntegerVector(x)) if (id < 1 || id > n) return false;
  return true;
}

bool linkage(SEXP x, bool anomer = false) {
  for (SEXP item : CharacterVector(x)) {
    const char* s = CHAR(item);
    int len = LENGTH(item);
    if (len < 2 || (s[0] != 'a' && s[0] != 'b' && s[0] != '?') ||
        (s[1] != '1' && s[1] != '2' && s[1] != '?')) return false;
    if (anomer) { if (len != 2) return false; continue; }
    if (len < 4 || s[2] != '-') return false;
    if (len == 4 && s[3] == '?') continue;
    if (len % 2 != 0) return false;
    for (int j = 3; j < len; ++j)
      if ((j % 2 && (s[j] < '1' || s[j] > '9')) || (!(j % 2) && s[j] != '/')) return false;
  }
  return true;
}

// Read only the checked v4 serialization schema, including graphs whose igraph
// C index has already been materialized. Extra/attributed fields use the public
// accessor path so their ordering and arbitrary user metadata are preserved.
struct GraphRecord {
  SEXP mono, sub, from, to, linkage, anomer, alditol, attributes;
  SEXP parts = R_NilValue, subs = R_NilValue;
};

bool transform_record(SEXP x, GraphRecord& a) {
  if (TYPEOF(x) != VECSXP || Rf_xlength(x) != 10 || !Rf_inherits(x, "igraph")) return false;
#if R_VERSION >= R_Version(4, 6, 0)
  if (R_getAttribCount(x) != 1) return false;
#else
  if (CDR(ATTRIB(x)) != R_NilValue) return false;
#endif
  SEXP cls = Rf_getAttrib(x, R_ClassSymbol);
  if (!strings(cls, 1) || std::strcmp(CHAR(STRING_ELT(cls, 0)), "igraph")) return false;
  List g(x);
  if (!plain(g[0], REALSXP, 1) || !plain(g[1], LGLSXP, 1) ||
      as<int>(g[1]) != TRUE || TYPEOF(g[9]) != ENVSXP) return false;
  double count = as<double>(g[0]);
  if (!R_FINITE(count) || count < 1 || count > 2048 || count != std::floor(count)) return false;
  int n = count;
  Environment env = g[9];
  if (!env.exists(".__igraph_version__.")) return false;
  SEXP version = env[".__igraph_version__."];
  if (!plain(version, INTSXP, 1) || INTEGER(version)[0] != 4) return false;
  for (int i = 4; i < 8; ++i) if (!Rf_isNull(g[i])) return false;
  if (TYPEOF(g[2]) != REALSXP || TYPEOF(g[3]) != REALSXP ||
      !no_attributes(g[2]) || !no_attributes(g[3]) ||
      Rf_xlength(g[2]) != Rf_xlength(g[3]) || Rf_xlength(g[2]) >= n) return false;
  if (TYPEOF(g[8]) != VECSXP || Rf_xlength(g[8]) != 4) return false;
  List attrs = g[8];
  if (!plain(attrs[0], REALSXP, 3)) return false;
  NumericVector flags = attrs[0];
  if (flags[0] != 1 || flags[1] != 0 || flags[2] != 1) return false;
  if (!fields(attrs[2], {"name", "mono", "sub"}) || !fields(attrs[3], {"linkage"})) return false;
  List vertices = attrs[2], edges = attrs[3];
  if (!strings(vertices[0], n) || !strings(vertices[1], n) || !strings(vertices[2], n) ||
      !strings(edges[0], Rf_xlength(g[2])) || !linkage(edges[0])) return false;
  CharacterVector names = vertices[0];
  for (int i = 0; i < n; ++i)
    if (as<std::string>(names[i]) != std::to_string(i + 1)) return false;
  SEXP ga = attrs[1];
  if (!fields(ga, {"anomer", "alditol"}) &&
      !fields(ga, {"anomer", "alditol", "floating_parts"}) &&
      !fields(ga, {"anomer", "alditol", "floating_substituents"}) &&
      !fields(ga, {"anomer", "alditol", "floating_parts", "floating_substituents"})) return false;
  List graph_attrs(ga);
  if (!strings(graph_attrs[0], 1) || !linkage(graph_attrs[0], true) ||
      !plain(graph_attrs[1], LGLSXP, 1) || as<int>(graph_attrs[1]) == NA_LOGICAL) return false;
  NumericVector from = g[2], to = g[3];
  for (int i = 0; i < from.size(); ++i) {
    if (!R_FINITE(from[i]) || !R_FINITE(to[i]) || from[i] < 0 || from[i] >= n ||
        to[i] < 0 || to[i] >= n || from[i] != std::floor(from[i]) || to[i] != std::floor(to[i]))
      return false;
  }
  a.mono = vertices[1]; a.sub = vertices[2];
  a.from = g[2]; a.to = g[3]; a.linkage = edges[0];
  a.anomer = graph_attrs[0]; a.alditol = graph_attrs[1];
  a.attributes = ga;
  if (graph_attrs.containsElementNamed("floating_parts")) {
    SEXP parts = graph_attrs["floating_parts"];
    if (TYPEOF(parts) != VECSXP || Rf_xlength(parts) == 0) return false;
    for (SEXP item : List(parts)) {
      if (!fields(item, {"root", "nodes", "linkage", "parents"})) return false;
      List p(item);
      if (!plain(p[0], INTSXP, 1) || !indices(p[0], n) ||
          !indices(p[1], n, true) || !strings(p[2], 1) || !linkage(p[2]) ||
          !indices(p[3], n)) return false;
    }
    a.parts = parts;
  }
  if (graph_attrs.containsElementNamed("floating_substituents")) {
    SEXP subs = graph_attrs["floating_substituents"];
    if (TYPEOF(subs) != VECSXP || Rf_xlength(subs) == 0) return false;
    for (SEXP item : List(subs)) {
      if (!fields(item, {"substituent", "parents"})) return false;
      List s(item);
      if (!strings(s[0], 1) || !indices(s[1], n)) return false;
    }
    a.subs = subs;
  }
  return true;
}

List forest_graph(const glyrepr_compact::Forest& f, List metadata, List model,
                  const GraphRecord& original, const std::string& operation,
                  std::unordered_map<int, CharacterVector>& node_names,
                  std::unordered_map<std::string, List>& plain_attributes,
                  CharacterVector vertex_fields, CharacterVector edge_fields) {
  const auto& t = f.tree;
  int n = t.mono.size(), ne = t.edge_nodes.size();
  bool same_vertices = true;
  for (int i = 0; i < n; ++i) if (t.vertices[i] != i) { same_vertices = false; break; }
  bool same_edges = same_vertices && ne == Rf_xlength(original.from) &&
    t.edge_nodes == f.edge_order;
  bool reuse_mono = same_vertices && operation != "generic";
  bool reuse_sub = same_vertices && operation != "remove_substituents" &&
    Rf_isNull(original.subs);
  bool reuse_linkage = same_edges && operation != "remove_linkages";
  // Only share ordinary endpoint vectors. A live igraph's ALTREP endpoints can
  // depend on its native index, so read those into independent plain vectors.
  bool reuse_endpoints = same_edges && !ALTREP(original.from) && !ALTREP(original.to);
  CharacterVector mono = reuse_mono ? CharacterVector(original.mono) : CharacterVector(n);
  CharacterVector sub = reuse_sub ? CharacterVector(original.sub) : CharacterVector(n);
  CharacterVector linkage = reuse_linkage ? CharacterVector(original.linkage) : CharacterVector(ne);
  NumericVector from = reuse_endpoints ? NumericVector(original.from) : NumericVector(ne);
  NumericVector to = reuse_endpoints ? NumericVector(original.to) : NumericVector(ne);
  std::vector<int> inverse(n);
  for (int i = 0; i < n; ++i) {
    int v = t.vertices[i];
    inverse[v] = i;
    if (!reuse_mono) mono[i] = t.mono[v];
    if (!reuse_sub) sub[i] = t.sub[v];
  }
  for (int i = 0; i < ne; ++i) {
    int v = t.edge_nodes[i];
    if (!reuse_endpoints) { from[i] = inverse[t.parent[v]]; to[i] = inverse[v]; }
    if (!reuse_linkage) linkage[i] = t.link[v];
  }
  auto entry = node_names.find(n);
  if (entry == node_names.end()) {
    CharacterVector names(n);
    for (int j = 0; j < n; ++j) names[j] = std::to_string(j+1);
    entry = node_names.emplace(n, names).first;
  }
  List graph_attrs;
  if (!metadata.size() && Rf_isNull(original.parts) && Rf_isNull(original.subs) &&
      operation != "remove_linkages") {
    graph_attrs = List(original.attributes);
  } else if (!metadata.size()) {
    // Immutable scalar metadata may be shared within this batch, just like
    // vertex names. igraph's attribute setters retain copy-on-modify behavior.
    std::string key = t.anomer + (t.alditol ? "1" : "0");
    auto found = plain_attributes.find(key);
    if (found == plain_attributes.end()) {
      found = plain_attributes.emplace(key,
        List::create(_["anomer"] = t.anomer, _["alditol"] = t.alditol)).first;
    }
    graph_attrs = found->second;
  } else {
    graph_attrs = List::create(_["anomer"] = t.anomer, _["alditol"] = t.alditol);
    if (metadata.containsElementNamed("floating_parts")) graph_attrs["floating_parts"] = metadata["floating_parts"];
    if (metadata.containsElementNamed("floating_substituents")) graph_attrs["floating_substituents"] = metadata["floating_substituents"];
  }
  List vertex_attrs = List::create(entry->second, mono, sub);
  vertex_attrs.attr("names") = vertex_fields;
  List edge_attrs = List::create(linkage);
  edge_attrs.attr("names") = edge_fields;
  List attrs = List::create(NumericVector::create(1, 0, 1), graph_attrs,
    vertex_attrs, edge_attrs);
  return graph_from_endpoints(from, to, n, true, attrs, model);
}
}

// [[Rcpp::export(name = ".compact_transform_graphs_native")]]
SEXP compact_transform_graphs_native(List graphs, CharacterVector keys,
    const std::string& operation, CharacterVector from, CharacterVector to,
    Function order, Function bliss, bool byte_order, SEXP model) {
  if (!supported_template(model) || graphs.size() != keys.size()) return R_NilValue;
  if (operation != "generic" && operation != "remove_linkages" &&
      operation != "remove_substituents") return R_NilValue;
  std::unordered_map<std::string, std::string> mapping;
  for (int i = 0; i < from.size(); ++i)
    mapping.emplace(as<std::string>(from[i]), as<std::string>(to[i]));
  List out(graphs.size());
  CharacterVector new_keys(keys.size());
  LogicalVector supported(graphs.size(), false);
  std::unordered_map<int, CharacterVector> node_names;
  std::unordered_map<std::string, List> plain_attributes;
  CharacterVector vertex_fields = CharacterVector::create("name", "mono", "sub");
  CharacterVector edge_fields = CharacterVector::create("linkage");
  std::unordered_map<std::string, int> seen;
  glyrepr_compact::LexicalCache lexical_cache;
  for (int i = 0; i < graphs.size(); ++i) {
    if (i % 128 == 0) checkUserInterrupt();
    try {
      GraphRecord a;
      if (!transform_record(graphs[i], a)) continue;
      bool changed = false;
      if (operation == "generic") {
        for (SEXP mono : CharacterVector(a.mono)) {
          auto found = mapping.find(CHAR(mono));
          if (found != mapping.end() && found->second != CHAR(mono)) { changed = true; break; }
        }
      } else if (operation == "remove_substituents") {
        for (SEXP sub : CharacterVector(a.sub)) if (LENGTH(sub)) { changed = true; break; }
        changed = changed || !Rf_isNull(a.subs);
      } else {
        changed = as<std::string>(a.anomer) != "??";
        for (SEXP link : CharacterVector(a.linkage))
          if (std::strcmp(CHAR(link), "?" "?-?")) { changed = true; break; }
        if (!Rf_isNull(a.parts))
          for (SEXP item : List(a.parts))
            if (as<std::string>(List(item)["linkage"]) != "?" "?-?") { changed = true; break; }
      }
      if (!changed) {
        new_keys[i] = keys[i]; out[i] = graphs[i]; supported[i] = true;
        seen.emplace(as<std::string>(keys[i]), i);
        continue;
      }
      auto f = glyrepr_compact::import_graph_data(a.mono, a.sub, a.anomer, a.alditol,
        a.from, a.linkage, a.parts, a.subs, false, a.to);
      f.tree.order_cache = &lexical_cache;
      glyrepr_compact::transform_graph(f, operation, mapping);
      // For connected trees with distinct sibling linkage ranks, neither
      // backbone selection nor branch order depends on subtree signatures.
      // Keep the full collation path for every tie and floating constraint.
      if (f.parts.empty() && f.subs.empty()) {
        f.tree.distinct_link_ranks = true;
        for (const auto& children : f.tree.children) {
          for (size_t j = 0; j < children.size(); ++j)
            for (size_t k = j + 1; k < children.size(); ++k)
              if (glyrepr_compact::rank_link(f.tree.link[children[j]]) ==
                  glyrepr_compact::rank_link(f.tree.link[children[k]]))
                f.tree.distinct_link_ranks = false;
          if (!f.tree.distinct_link_ranks) break;
        }
      }
      if (!f.parts.empty() || !f.subs.empty()) glyrepr_compact::validate_forest(f);
      List result = glyrepr_compact::finish_forest(f, order, byte_order, bliss, false);
      const std::string& key = f.tree.sequence;
      new_keys[i] = key;
      auto found = seen.find(key);
      if (found == seen.end()) {
        out[i] = forest_graph(f, result, model, a, operation, node_names, plain_attributes,
                              vertex_fields, edge_fields);
        seen.emplace(key, i);
      } else out[i] = out[found->second];
      supported[i] = true;
    } catch (const std::exception&) {
      // The R reference path retains diagnostics and unusual graph semantics.
    }
  }
  return List::create(_["keys"] = new_keys, _["graphs"] = out, _["supported"] = supported);
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
