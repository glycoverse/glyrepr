#include <Rcpp.h>
#include <algorithm>
#include <string>
#include <vector>
using namespace Rcpp;

// Validated ordinary trees only. R owns graph extraction, metadata and errors.
// Retain R collation/stable tie handling through base sort/order at branch points.
// [[Rcpp::export]]
List sequence_cpp(IntegerMatrix endpoints, CharacterVector mono,
                  CharacterVector sub, CharacterVector linkage) {
  int n = mono.size(), m = endpoints.nrow();
  if (n < 1 || endpoints.ncol() != 2 || sub.size() != n ||
      linkage.size() != m || m != n - 1) stop("Expected a nonempty tree");
  std::vector<std::vector<int>> children(n), ordered(n);
  std::vector<int> indegree(n, 0), target(m), depths(n, 0), traversal;
  std::vector<std::string> labels(n), links(m), signatures(n);
  std::vector<double> ranks(m, 1.0);
  for (int v = 0; v < n; ++v) {
    labels[v] = as<std::string>(mono[v]);
    if (sub[v] != NA_STRING) {
      std::string s = as<std::string>(sub[v]);
      s.erase(std::remove(s.begin(), s.end(), ','), s.end());
      labels[v] += s;
    }
  }
  for (int e = 0; e < m; ++e) {
    int p = endpoints(e, 0) - 1, c = endpoints(e, 1) - 1;
    if (p < 0 || p >= n || c < 0 || c >= n || ++indegree[c] > 1)
      stop("Invalid tree endpoints");
    children[p].push_back(e); target[e] = c;
    links[e] = as<std::string>(linkage[e]);
    std::string pos = links[e].substr(links[e].rfind('-') + 1);
    if (pos != "?" && pos.find('/') == std::string::npos)
      ranks[e] = 1.0 / std::stod(pos);
  }
  int root = -1;
  for (int v = 0; v < n; ++v) if (!indegree[v]) {
    if (root != -1) stop("Multiple roots");
    root = v;
  }
  if (root < 0) stop("Missing root");
  traversal.push_back(root);
  for (size_t i = 0; i < traversal.size(); ++i) {
    if (i % 1024 == 0) checkUserInterrupt();
    for (int e : children[traversal[i]]) traversal.push_back(target[e]);
  }
  if (static_cast<int>(traversal.size()) != n) stop("Disconnected tree");
  Environment base = Environment::base_env();
  Function rsort = base["sort"], rorder = base["order"];
  for (auto it = traversal.rbegin(); it != traversal.rend(); ++it) {
    int v = *it, k = children[v].size();
    signatures[v] = labels[v];
    if (!k) continue;
    CharacterVector tokens(k), child_sigs(k);
    IntegerVector child_depths(k);
    NumericVector child_ranks(k);
    for (int j = 0; j < k; ++j) {
      int e = children[v][j], c = target[e];
      depths[v] = std::max(depths[v], depths[c] + 1);
      tokens[j] = links[e] + "->" + signatures[c];
      child_sigs[j] = signatures[c];
      child_depths[j] = depths[c]; child_ranks[j] = ranks[e];
    }
    CharacterVector sorted = k > 1 ? as<CharacterVector>(rsort(tokens)) : tokens;
    signatures[v] += "{";
    for (int j = 0; j < k; ++j) {
      if (j) signatures[v] += ",";
      signatures[v] += as<std::string>(sorted[j]);
    }
    signatures[v] += "}";
    if (k == 1) { ordered[v] = children[v]; continue; }
    IntegerVector back = rorder(child_depths, child_ranks, child_sigs,
                                Named("decreasing") = true);
    IntegerVector branches = rorder(child_ranks, child_sigs,
                                    Named("decreasing") = true);
    int backbone = back[0] - 1;
    ordered[v].push_back(children[v][backbone]);
    for (int j : branches) if (j - 1 != backbone)
      ordered[v].push_back(children[v][j - 1]);
  }
  // Explicit event stack avoids recursive vector/string concatenation.
  struct Event { int kind; int id; bool branch; };
  std::vector<Event> stack{{0, root, false}};
  NumericVector vertices(n), edges(m);
  int vi = 0, ei = 0;
  std::string out;
  while (!stack.empty()) {
    Event event = stack.back(); stack.pop_back();
    if (event.kind == 1) {
      vertices[vi++] = event.id + 1; out += labels[event.id];
    } else if (event.kind == 2) {
      edges[ei++] = event.id + 1; out += "(" + links[event.id] + ")";
      if (event.branch) out += "]";
    } else if (event.kind == 3) {
      out += "[";
    } else {
      int v = event.id;
      stack.push_back({1, v, false});
      for (int j = static_cast<int>(ordered[v].size()) - 1; j >= 0; --j) {
        int e = ordered[v][j];
        stack.push_back({2, e, j > 0});
        stack.push_back({0, target[e], false});
        if (j > 0) stack.push_back({3, 0, false});
      }
    }
  }
  return List::create(Named("vertices") = vertices, Named("edges") = edges,
                      Named("iupac") = out);
}
