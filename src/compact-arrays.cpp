#include "compact-core.h"

namespace glyrepr_compact {

// R validates field types, lengths and chemistry before crossing this boundary.
// This importer verifies topology without constructing an intermediate igraph.
Forest import_arrays(const List& a) {
  Forest f;
  Tree& t = f.tree;
  t.mono = as<std::vector<std::string>>(a["mono"]);
  t.sub = as<std::vector<std::string>>(a["sub"]);
  t.anomer = as<std::string>(a["anomer"]);
  t.alditol = as<bool>(a["alditol"]);
  int n = t.mono.size();
  t.parent.assign(n, -1);
  t.link.resize(n);
  t.children.resize(n);
  IntegerVector edges = a["edges"];
  CharacterVector links = a["linkage"];
  for (int i = 0; i < links.size(); ++i) {
    int parent = edges[2*i]-1, child = edges[2*i+1]-1;
    if (parent == child || t.parent[child] != -1)
      throw std::runtime_error("self edge or multiple parents");
    t.parent[child] = parent;
    t.children[parent].push_back(child);
    t.link[child] = as<std::string>(links[i]);
    f.edge_order.push_back(child);
  }
  // Iterative traversal also rejects rootless cycles before recursive emit/cache.
  std::vector<int> component(n, -1), roots;
  for (int i = 0; i < n; ++i) if (t.parent[i] == -1) roots.push_back(i);
  int visited = 0;
  for (int root : roots) {
    std::vector<int> pending = {root};
    while (!pending.empty()) {
      int v = pending.back(); pending.pop_back();
      if (component[v] != -1) throw std::runtime_error("cyclic graph");
      component[v] = root; ++visited;
      for (int child : t.children[v]) pending.push_back(child);
    }
  }
  if (visited != n) throw std::runtime_error("glycan components must be out trees");
  auto indices = [](SEXP x) {
    auto ids = as<std::vector<int>>(x);
    for (int& id : ids) --id;
    std::sort(ids.begin(), ids.end());
    return ids;
  };
  std::set<int> declared;
  if (a.containsElementNamed("floating_parts") && !Rf_isNull(a["floating_parts"])) {
    List parts = a["floating_parts"];
    for (SEXP item : parts) {
      List p(item);
      Part part;
      part.root = as<int>(p["root"])-1;
      part.nodes = indices(p["nodes"]);
      part.parents = indices(p["parents"]);
      part.linkage = as<std::string>(p["linkage"]);
      if (t.parent[part.root] != -1 || !declared.insert(part.root).second)
        throw std::runtime_error("floating roots must identify distinct components");
      std::vector<int> actual;
      for (int i = 0; i < n; ++i) if (component[i] == part.root) actual.push_back(i);
      if (actual != part.nodes) throw std::runtime_error("floating nodes do not match component");
      for (int parent : part.parents) if (component[parent] == part.root)
        throw std::runtime_error("floating parent cannot belong to its own component");
      f.parts.push_back(part);
    }
  }
  for (int root : roots) if (!declared.count(root)) {
    if (f.root != -1) throw std::runtime_error("undeclared disconnected component");
    f.root = root;
  }
  if (f.root == -1) throw std::runtime_error("missing main component");
  if (a.containsElementNamed("floating_substituents") && !Rf_isNull(a["floating_substituents"])) {
    List subs = a["floating_substituents"];
    for (SEXP item : subs) {
      List s(item);
      f.subs.push_back({as<std::string>(s["substituent"]), indices(s["parents"])});
    }
  }
  return f;
}

} // namespace glyrepr_compact

// [[Rcpp::export(name = ".compact_arrays_native")]]
List compact_arrays_native(List records, Function order, Function bliss, bool byte_order = false) {
  List out(records.size());
  for (int i = 0; i < records.size(); ++i) {
    if (i % 128 == 0) checkUserInterrupt();
    if (Rf_isNull(records[i])) {out[i] = List::create(_["status"]="missing"); continue;}
    try {
      List a = records[i];
      CharacterVector mono = a["mono"];
      if (mono.size() > 2048) {
        out[i] = List::create(_["status"]="unsupported", _["reason"]="native size guard");
        continue;
      }
      auto f = glyrepr_compact::import_arrays(a);
      if (!f.parts.empty() || !f.subs.empty()) glyrepr_compact::validate_forest(f);
      out[i] = glyrepr_compact::finish_forest(f, order, byte_order, bliss);
    } catch (const std::exception& e) {
      out[i] = List::create(_["status"]="invalid", _["reason"]=e.what());
    }
  }
  return out;
}
