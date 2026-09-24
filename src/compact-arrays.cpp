#include "compact-core.h"

namespace glyrepr_compact {

// R validates field types, lengths and chemistry before crossing this boundary.
// This importer verifies topology without constructing an intermediate igraph.
Forest import_arrays(const List& a, bool preserve_parent_order = false) {
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
  auto indices = [](SEXP x, bool sorted = true) {
    auto ids = as<std::vector<int>>(x);
    for (int& id : ids) --id;
    if (sorted) std::sort(ids.begin(), ids.end());
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
      part.parents = indices(p["parents"], !preserve_parent_order);
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
      f.subs.push_back({as<std::string>(s["substituent"]), indices(s["parents"], !preserve_parent_order)});
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


namespace glyrepr_compact {

// Serialization preserves the supplied metadata IDs and constraint order.
// It deliberately does not resolve singleton domains or canonicalize a graph.
std::string graph_sequence(Forest& f, Function& order, bool byte_order) {
  Tree& t = f.tree;
  cache_tree(t, order, byte_order);
  std::string result;
  for (const auto& s : f.subs) {
    auto parents = s.parents;
    for (int& id : parents) ++id;
    result += "{" + s.token + (parents.empty() ? "" : "|" + int_join(parents)) + "}";
  }
  for (const auto& p : f.parts) {
    auto parents = p.parents;
    for (int& id : parents) ++id;
    result += "{" + emit(t, p.root, order, byte_order) + "(" + p.linkage + ")" +
      (parents.empty() ? "" : "|" + int_join(parents)) + "}";
  }
  return result + emit(t, f.root, order, byte_order) +
    (t.alditol ? "-ol" : "") + "(" + t.anomer + "-";
}

List finish_graph(Forest& f, Function& order, bool byte_order, Function& bliss) {
  std::vector<int> original_edges(f.tree.mono.size(), NA_INTEGER);
  for (int i = 0; i < (int)f.edge_order.size(); ++i)
    original_edges[f.edge_order[i]] = i + 1;
  List out = finish_forest(f, order, byte_order, bliss);
  IntegerVector vertices(f.tree.vertices.size()), edges(f.tree.edge_nodes.size());
  for (int i = 0; i < vertices.size(); ++i) vertices[i] = f.tree.vertices[i] + 1;
  for (int i = 0; i < edges.size(); ++i) edges[i] = original_edges[f.tree.edge_nodes[i]];
  out["vertex_order"] = vertices;
  out["edge_order"] = edges;
  return out;
}

void transform_graph(Forest& f, const std::string& operation,
                     const std::unordered_map<std::string, std::string>& mapping) {
  Tree& t = f.tree;
  if (operation == "generic") {
    for (auto& mono : t.mono) {
      auto found = mapping.find(mono);
      if (found != mapping.end()) mono = found->second;
    }
  } else if (operation == "remove_linkages") {
    for (int node : f.edge_order) t.link[node] = "??-?";
    t.anomer = "??";
    for (auto& part : f.parts) part.linkage = "??-?";
  } else if (operation == "remove_substituents") {
    std::fill(t.sub.begin(), t.sub.end(), "");
    f.subs.clear();
  } else if (operation == "fill_anomer_pos") {
    auto fill = [&](std::string& value, int node) {
      if (value.size() > 1 && value[1] == '?') {
        auto found = mapping.find(t.mono[node]);
        if (found != mapping.end()) value.replace(1, 1, found->second);
      }
    };
    fill(t.anomer, f.root);
    for (int node : f.edge_order) fill(t.link[node], node);
    for (auto& part : f.parts) fill(part.linkage, part.root);
  }
}

} // namespace glyrepr_compact

// Records have checked scalar types and indices before crossing this boundary.
// Errors are returned to R, which replays the public reference path for diagnostics.
// [[Rcpp::export(name = ".compact_graphs_native")]]
List compact_graphs_native(List records, std::string mode, bool validate,
                           std::string operation, CharacterVector from,
                           CharacterVector to, Function order, Function bliss,
                           bool byte_order = false) {
  std::unordered_map<std::string, std::string> mapping;
  for (int i = 0; i < from.size(); ++i)
    mapping.emplace(as<std::string>(from[i]), as<std::string>(to[i]));
  List out(records.size());
  for (int i = 0; i < records.size(); ++i) {
    if (i % 128 == 0) checkUserInterrupt();
    if (Rf_isNull(records[i])) continue;
    try {
      List a = records[i];
      CharacterVector mono = a["mono"];
      if (mono.size() > 2048) throw glyrepr_compact::Unsupported("native size guard");
      auto f = glyrepr_compact::import_arrays(a, mode == "serialize");
      if (validate) glyrepr_compact::validate_forest(f);
      if (mode == "validate") {
        out[i] = List::create(_["status"] = "ok");
      } else if (mode == "serialize") {
        out[i] = List::create(_["status"] = "ok",
          _["iupac"] = glyrepr_compact::graph_sequence(f, order, byte_order));
      } else {
        glyrepr_compact::transform_graph(f, operation, mapping);
        if (!operation.empty() && (!f.parts.empty() || !f.subs.empty()))
          glyrepr_compact::validate_forest(f);
        out[i] = glyrepr_compact::finish_graph(f, order, byte_order, bliss);
      }
    } catch (const std::exception& e) {
      out[i] = List::create(_["status"] = "fallback", _["reason"] = e.what());
    }
  }
  return out;
}

// [[Rcpp::export(name = ".compact_localizations_native")]]
List compact_localizations_native(List record, List domains, int combinations) {
  try {
    auto original = glyrepr_compact::import_arrays(record);
    if (original.tree.mono.size() > 2048)
      return List::create(_["status"] = "fallback");
    int np = original.parts.size(), ns = original.subs.size();
    if (domains.size() != np + ns) stop("invalid domain count");
    std::vector<std::vector<int>> candidates;
    for (SEXP domain : domains) candidates.push_back(as<std::vector<int>>(domain));
    List records, assignments;
    for (int combination = 0; combination < combinations; ++combination) {
      if (combination % 128 == 0) checkUserInterrupt();
      auto f = original;
      int remainder = combination;
      IntegerVector selected(domains.size());
      for (int i = 0; i < domains.size(); ++i) {
        if (candidates[i].empty()) stop("empty candidate domain");
        int parent = candidates[i][remainder % candidates[i].size()];
        remainder /= candidates[i].size();
        selected[i] = parent;
        if (i < np) f.parts[i].parents = {parent - 1};
        else f.subs[i - np].parents = {parent - 1};
      }
      try {
        glyrepr_compact::validate_forest(f);
      } catch (const std::exception&) {
        continue;
      }
      // Materialize accepted assignments in original vertex and edge order.
      // No singleton normalization or pruning is necessary for complete domains.
      for (const auto& p : f.parts) {
        int parent = p.parents[0];
        f.tree.parent[p.root] = parent;
        f.tree.children[parent].push_back(p.root);
        f.tree.link[p.root] = p.linkage;
        f.edge_order.push_back(p.root);
      }
      for (const auto& s : f.subs) {
        int parent = s.parents[0];
        auto tokens = glyrepr_compact::split(f.tree.sub[parent], ',');
        tokens.push_back(s.token);
        f.tree.sub[parent] = glyrepr_compact::collapse_subs(tokens);
      }
      f.tree.vertices.resize(f.tree.mono.size());
      std::iota(f.tree.vertices.begin(), f.tree.vertices.end(), 0);
      f.tree.edge_nodes = f.edge_order;
      records.push_back(glyrepr_compact::serialize_tree(f.tree));
      assignments.push_back(selected);
    }
    return List::create(_["status"] = "ok", _["records"] = records,
                        _["assignments"] = assignments);
  } catch (const std::exception& e) {
    return List::create(_["status"] = "fallback", _["reason"] = e.what());
  }
}

// [[Rcpp::export(name = ".compact_localize_parts_native")]]
List compact_localize_parts_native(List record, IntegerVector part_ids,
                                   IntegerVector parents) {
  try {
    auto f = glyrepr_compact::import_arrays(record, true);
    if (f.tree.mono.size() > 2048) return List::create(_["status"] = "fallback");
    auto selected = f;
    std::vector<bool> attached(f.parts.size(), false);
    for (int i = 0; i < part_ids.size(); ++i) {
      int id = part_ids[i] - 1;
      if (id < 0 || id >= (int)f.parts.size() || parents[i] < 1 ||
          parents[i] > (int)f.tree.mono.size()) stop("invalid assignment");
      selected.parts[id].parents = {parents[i] - 1};
      attached[id] = true;
    }
    glyrepr_compact::validate_forest(selected);
    for (int i = 0; i < part_ids.size(); ++i) {
      const auto& p = f.parts[part_ids[i] - 1];
      int parent = parents[i] - 1;
      f.tree.parent[p.root] = parent;
      f.tree.children[parent].push_back(p.root);
      f.tree.link[p.root] = p.linkage;
      f.edge_order.push_back(p.root);
    }
    std::vector<glyrepr_compact::Part> remaining;
    for (int i = 0; i < (int)f.parts.size(); ++i) if (!attached[i]) {
      auto p = f.parts[i];
      p.nodes.clear();
      std::vector<int> pending = {p.root};
      while (!pending.empty()) {
        int node = pending.back(); pending.pop_back();
        p.nodes.push_back(node);
        for (int child : f.tree.children[node]) pending.push_back(child);
      }
      std::sort(p.nodes.begin(), p.nodes.end());
      if (!p.parents.empty()) {
        std::vector<int> keep;
        for (int parent : p.parents)
          if (!std::binary_search(p.nodes.begin(), p.nodes.end(), parent)) keep.push_back(parent);
        if (keep.empty()) stop("empty remaining domain");
        p.parents = keep;
      }
      remaining.push_back(p);
    }
    f.parts = remaining;
    glyrepr_compact::validate_forest(f);
    f.tree.vertices.resize(f.tree.mono.size());
    std::iota(f.tree.vertices.begin(), f.tree.vertices.end(), 0);
    f.tree.edge_nodes = f.edge_order;
    List out = glyrepr_compact::serialize_tree(f.tree);
    List parts(f.parts.size());
    for (int i = 0; i < (int)f.parts.size(); ++i) {
      const auto& p = f.parts[i];
      auto nodes = p.nodes, domains = p.parents;
      for (int& id : nodes) ++id;
      for (int& id : domains) ++id;
      parts[i] = List::create(_["root"] = p.root + 1, _["nodes"] = wrap(nodes),
                             _["linkage"] = p.linkage, _["parents"] = wrap(domains));
    }
    out["floating_parts"] = parts;
    return out;
  } catch (const std::exception& e) {
    return List::create(_["status"] = "fallback", _["reason"] = e.what());
  }
}
