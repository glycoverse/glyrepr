#include <Rcpp.h>
#include <algorithm>
#include <regex>
#include <unordered_map>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

// No igraph API or igraph object layout is used here. Indices stay in arrays.
struct Tree {
  std::vector<std::string> mono, link, sig;
  std::vector<int> parent, depth, vertices, edge_nodes;
  std::vector<std::vector<int>> children;
  std::string anomer, sequence;
  bool alditol = false;
};

std::vector<int> lexical_order(const std::vector<std::string>& strings,
                               Function& order, bool byte_order) {
  std::vector<int> ids(strings.size());
  std::iota(ids.begin(), ids.end(), 0);
  if (byte_order) {
    std::stable_sort(ids.begin(), ids.end(), [&](int a, int b) {
      return strings[a] < strings[b];
    });
  } else if (strings.size() > 1) {
    // Preserve R's active collation, rather than silently imposing ASCII order.
    IntegerVector p = order(wrap(strings));
    for (int i = 0; i < p.size(); ++i) ids[i] = p[i] - 1;
  }
  return ids;
}

double rank_link(const std::string& link) {
  auto p = link.substr(3);
  return p == "?" || p.find('/') != std::string::npos ? 1.0 : 1.0 / (p[0] - '0');
}

std::string emit(Tree& t, int node, Function& order, bool byte_order) {
  auto kids = t.children[node];
  if (!kids.empty()) {
    std::vector<std::string> signatures;
    for (int k : kids) signatures.push_back(t.sig[k]);
    auto lexical = lexical_order(signatures, order, byte_order);
    std::vector<int> ranks(kids.size());
    int r = 0;
    for (int j = 0; j < (int)lexical.size(); ++j) {
      if (j && signatures[lexical[j]] != signatures[lexical[j-1]]) ++r;
      ranks[lexical[j]] = r;
    }
    std::vector<int> ids(kids.size());
    std::iota(ids.begin(), ids.end(), 0);
    auto better = [&](int a, int b) {
      double ra = rank_link(t.link[kids[a]]), rb = rank_link(t.link[kids[b]]);
      if (ra != rb) return ra > rb;
      return ranks[a] > ranks[b];
    };
    int backbone = 0;
    for (int j = 1; j < (int)kids.size(); ++j) {
      if (t.depth[kids[j]] > t.depth[kids[backbone]] ||
          (t.depth[kids[j]] == t.depth[kids[backbone]] && better(j, backbone))) backbone = j;
    }
    std::stable_sort(ids.begin(), ids.end(), better);
    std::string result = emit(t, kids[backbone], order, byte_order) + "(" + t.link[kids[backbone]] + ")";
    t.edge_nodes.push_back(kids[backbone]);
    for (int j : ids) if (j != backbone) {
      result += "[" + emit(t, kids[j], order, byte_order) + "(" + t.link[kids[j]] + ")]";
      t.edge_nodes.push_back(kids[j]);
    }
    t.vertices.push_back(node);
    return result + t.mono[node];
  }
  t.vertices.push_back(node);
  return t.mono[node];
}

List parse_tree(std::string s, const std::unordered_map<std::string, int>& known,
                Function& order, bool byte_order) {
  if (s.find_first_of("{}") != std::string::npos)
    return List::create(_["status"]="fallback", _["reason"]="floating metadata");
  if (s.empty() || std::any_of(s.begin(), s.end(), [](unsigned char c){ return std::isspace(c); }))
    throw std::runtime_error("empty input or whitespace");
  Tree t;
  static const std::regex reducing("\\(([ab?][12?])-$");
  static const std::regex linkage("^[ab?][12?]-([1-9](/[1-9])*|\\?)$");
  static const std::regex mono_pattern("^([DL]-)?([A-Za-z]|[0-9][A-Za-z])[A-Za-z0-9?/]*$");
  std::smatch match;
  if (std::regex_search(s, match, reducing)) {
    t.anomer = match[1]; s.resize(match.position());
  }
  if (s.size() >= 3 && s.substr(s.size()-3) == "-ol") {
    t.alditol = true; s.resize(s.size()-3);
  }
  if (s.find("-ol") != std::string::npos) throw std::runtime_error("misplaced alditol");
  // Read right to left, so every parent is already allocated.
  int cursor = (int)s.size()-1, current = -1;
  std::vector<int> stack;
  std::vector<bool> branch_has_node;
  while (cursor >= 0) {
    if (s[cursor] == ']') {
      if (current < 0) throw std::runtime_error("branch without parent");
      stack.push_back(current); branch_has_node.push_back(false); --cursor; continue;
    }
    if (s[cursor] == '[') {
      if (stack.empty() || !branch_has_node.back()) throw std::runtime_error("malformed branch");
      current = stack.back(); stack.pop_back(); branch_has_node.pop_back(); --cursor; continue;
    }
    std::string link;
    if (s[cursor] == ')') {
      int end = cursor--;
      while (cursor >= 0 && s[cursor] != '(') --cursor;
      if (cursor < 0) throw std::runtime_error("missing opening parenthesis");
      link = s.substr(cursor+1, end-cursor-1); --cursor;
      if (link == "?-?") link = "??-?";
      // Normalize any ambiguous acceptor containing '?' to unknown, like R.
      if (link.size() >= 4 && link[2] == '-' && link.substr(3).find('?') != std::string::npos) {
        static const std::regex unknown("^[ab?][12?]-([1-9]/)*\\?(/[1-9])*$");
        if (std::regex_match(link, unknown)) link = link.substr(0,3) + "?";
      }
      if (!std::regex_match(link, linkage)) throw std::runtime_error("invalid linkage");
    }
    int end = cursor;
    while (cursor >= 0 && s[cursor] != '[' && s[cursor] != ']' && s[cursor] != ')') --cursor;
    std::string mono = s.substr(cursor+1, end-cursor);
    if (!std::regex_match(mono, mono_pattern)) throw std::runtime_error("invalid residue syntax");
    auto it = known.find(mono);
    if (it == known.end())
      return List::create(_["status"]="fallback", _["reason"]="substituted or unknown residue");
    if ((current == -1) != link.empty()) throw std::runtime_error("missing or misplaced linkage");
    int id = t.mono.size();
    if (id >= 2048) return List::create(_["status"]="fallback", _["reason"]="size guard");
    t.mono.push_back(mono); t.link.push_back(link); t.parent.push_back(current);
    t.children.emplace_back(); t.depth.push_back(0); t.sig.emplace_back();
    if (current >= 0) t.children[current].push_back(id);
    else if (t.anomer.empty()) t.anomer = "?" + std::to_string(it->second);
    current = id;
    if (!branch_has_node.empty()) branch_has_node.back() = true;
  }
  if (!stack.empty() || t.mono.empty()) throw std::runtime_error("unbalanced branch or empty tree");
  for (int n = (int)t.mono.size()-1; n >= 0; --n) {
    int occupied = 0;
    std::vector<std::string> tokens;
    for (int k : t.children[n]) {
      auto pos = t.link[k].substr(3);
      if (pos.size() == 1 && pos != "?") {
        int bit = 1 << (pos[0]-'0');
        if (occupied & bit) throw std::runtime_error("duplicated linkage positions");
        occupied |= bit;
      }
      t.depth[n] = std::max(t.depth[n], t.depth[k]+1);
      tokens.push_back(t.link[k] + "->" + t.sig[k]);
    }
    t.sig[n] = t.mono[n];
    if (!tokens.empty()) {
      auto ids = lexical_order(tokens, order, byte_order);
      t.sig[n] += "{";
      for (int j = 0; j < (int)ids.size(); ++j) t.sig[n] += (j ? "," : "") + tokens[ids[j]];
      t.sig[n] += "}";
    }
  }
  t.sequence = emit(t, 0, order, byte_order) + (t.alditol ? "-ol" : "") + "(" + t.anomer + "-";
  int n = t.mono.size();
  CharacterVector monos(n), links(n-1);
  IntegerVector parents(n), edges(2*(n-1));
  std::vector<int> inverse(n);
  for (int j = 0; j < n; ++j) { inverse[t.vertices[j]] = j+1; monos[j] = t.mono[t.vertices[j]]; }
  for (int j = 0; j < n; ++j) {
    int p = t.parent[t.vertices[j]]; parents[j] = p < 0 ? 0 : inverse[p];
  }
  for (int j = 0; j < n-1; ++j) {
    int k = t.edge_nodes[j]; edges[2*j] = inverse[t.parent[k]]; edges[2*j+1] = inverse[k]; links[j] = t.link[k];
  }
  return List::create(_["status"]="ok", _["iupac"]=t.sequence, _["mono"]=monos,
                      _["parent"]=parents, _["edges"]=edges, _["linkage"]=links,
                      _["anomer"]=t.anomer, _["alditol"]=t.alditol);
}

// [[Rcpp::export]]
List compact_parse(CharacterVector strings, CharacterVector residues,
                   IntegerVector anomer_positions, Function order, bool byte_order=false) {
  std::unordered_map<std::string,int> known;
  for (int i = 0; i < residues.size(); ++i) known.emplace(as<std::string>(residues[i]), anomer_positions[i]);
  List out(strings.size());
  for (int i = 0; i < strings.size(); ++i) {
    if (i % 128 == 0) checkUserInterrupt();
    if (strings[i] == NA_STRING) { out[i] = List::create(_["status"]="missing"); continue; }
    try { out[i] = parse_tree(as<std::string>(strings[i]), known, order, byte_order); }
    catch (const std::exception& e) {
      out[i] = List::create(_["status"]="error", _["reason"]=e.what());
    }
  }
  return out;
}
