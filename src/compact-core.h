#ifndef GLYREPR_COMPACT_CORE_H
#define GLYREPR_COMPACT_CORE_H

#include <Rcpp.h>
#include <algorithm>
#include <cctype>
#include <climits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <regex>
#include <unordered_map>
#include <unordered_set>
#include <numeric>
#include <functional>
#include <set>
using namespace Rcpp;
namespace glyrepr_compact {

struct Unsupported : std::runtime_error {
  using std::runtime_error::runtime_error;
};

// Tree operations use arrays. Unresolved metadata uses one public BLISS callback;
// no igraph object layout is read and the final glycan graph is built afterward.
struct Tree {
  std::vector<std::string> mono, sub, link, sig;
  std::vector<int> parent, depth, vertices, edge_nodes, symmetry;
  std::vector<std::vector<int>> children;
  std::string anomer, sequence;
  bool alditol = false;
};


inline std::string residue_label(const Tree& t, int n) {
  std::string sub = t.sub[n];
  sub.erase(std::remove(sub.begin(), sub.end(), ','), sub.end());
  return t.mono[n] + sub;
}

inline std::vector<std::string> split(const std::string& s, char sep) {
  std::vector<std::string> out;
  size_t begin = 0;
  if (s.empty()) return out;
  for (size_t p = 0; p <= s.size(); ++p) if (p == s.size() || s[p] == sep) {
    out.push_back(s.substr(begin, p-begin)); begin = p+1;
  }
  return out;
}
inline std::string join(const std::vector<std::string>& x, const std::string& sep) {
  std::string out;
  for (size_t i = 0; i < x.size(); ++i) out += (i ? sep : "") + x[i];
  return out;
}
inline std::string int_join(const std::vector<int>& x) {
  std::vector<std::string> strings;
  for (int i : x) strings.push_back(std::to_string(i));
  return join(strings, ",");
}
const std::regex sub_regex("([1-9](/[1-9])*|\\?)(PPEtn|PEtn|NAc|NGc|Pyr|Me|Ac|PC|Gc|P|S|N)");
inline std::vector<int> positions(const std::string& s) {
  if (s.empty() || s[0] == '?') return {};
  std::vector<int> out;
  for (char c : s) { if (c >= '1' && c <= '9') out.push_back(c-'0'); else if (c != '/') break; }
  std::sort(out.begin(), out.end());
  out.erase(std::unique(out.begin(), out.end()), out.end());
  return out;
}
inline std::string normalize_sub(std::string s) {
  if (!std::regex_match(s, sub_regex)) throw std::runtime_error("invalid substituent");
  if (s[0] == '?') return s;
  size_t k = s.find_first_not_of("123456789/");
  auto pos = positions(s); std::vector<std::string> p;
  for (int i : pos) p.push_back(std::to_string(i));
  return join(p, "/") + s.substr(k);
}
inline bool matching(std::vector<std::vector<int>> domains) {
  std::sort(domains.begin(), domains.end(), [](const std::vector<int>& a, const std::vector<int>& b){return a.size() < b.size();});
  // Augmenting-path bipartite matching avoids exponential carbon-slot search.
  std::unordered_map<int,int> owners;
  std::function<bool(int,std::unordered_set<int>&)> assign = [&](int i, std::unordered_set<int>& seen) {
    for (int slot : domains[i]) if (seen.insert(slot).second) {
      auto it = owners.find(slot);
      if (it == owners.end() || assign(it->second, seen)) {owners[slot]=i; return true;}
    }
    return false;
  };
  for (int i=0; i<(int)domains.size(); ++i) {std::unordered_set<int> seen; if (!assign(i,seen)) return false;}
  return true;
}
inline std::string collapse_subs(std::vector<std::string> tokens) {
  std::vector<std::vector<int>> domains;
  for (auto& token : tokens) {
    token = normalize_sub(token);
    auto p = positions(token); if (!p.empty()) domains.push_back(p);
  }
  if (!matching(domains)) throw std::runtime_error("conflicting substituent positions");
  std::stable_sort(tokens.begin(), tokens.end(), [](const std::string& a, const std::string& b) {
    return (a[0]=='?' ? 99 : a[0]-'0') < (b[0]=='?' ? 99 : b[0]-'0');
  });
  return join(tokens, ",");
}
inline std::pair<std::string,std::string> extract_basic(std::string mono) {
  bool neu = mono.compare(0,3,"Neu") == 0;
  auto marker_pos = [&](const std::string& marker) {
    size_t p = 0;
    while ((p = mono.find(marker,p)) != std::string::npos) {
      if (p==0 || (mono[p-1]!='/' && !std::isdigit((unsigned char)mono[p-1]))) return p;
      ++p;
    }
    return std::string::npos;
  };
  size_t ac = marker_pos("5Ac"), gc = marker_pos("5Gc");
  if (neu && ac!=std::string::npos && gc!=std::string::npos) throw std::runtime_error("conflicting Neu markers");
  std::string base;
  if (neu && (ac!=std::string::npos || gc!=std::string::npos)) {
    base = mono.compare(0,4,"Neuf")==0 ? "Neuf" : "Neu";
    base += ac!=std::string::npos ? "5Ac" : "5Gc";
    mono.erase(ac!=std::string::npos ? ac : gc, 3);
  }
  std::vector<std::string> subs;
  for (std::sregex_iterator it(mono.begin(),mono.end(),sub_regex), end; it!=end; ++it) subs.push_back(it->str());
  std::string cleaned = mono;
  for (const auto& token : subs) {
    size_t p = cleaned.find(token); if (p!=std::string::npos) cleaned.erase(p,token.size());
  }
  return {base.empty() ? cleaned : base, collapse_subs(subs)};
}
inline std::pair<std::string,std::string> extract_residue(const std::string& mono,
    const std::unordered_map<std::string,int>& known,
    const std::unordered_map<std::string,std::string>& configurations) {
  if (known.count(mono)) return {mono,""};
  auto result = extract_basic(mono);
  if (!known.count(result.first) && mono.size()>2 && (mono[0]=='D'||mono[0]=='L') && mono[1]=='-') {
    auto unconfigured = extract_basic(mono.substr(2));
    auto it = configurations.find(unconfigured.first);
    if (it!=configurations.end() && it->second.substr(0,2)==mono.substr(0,2)) result={it->second,unconfigured.second};
  }
  if (!known.count(result.first)) throw std::runtime_error("unknown monosaccharide");
  return result;
}

inline std::vector<int> lexical_order(const std::vector<std::string>& strings,
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

inline double rank_link(const std::string& link) {
  auto p = link.substr(3);
  return p == "?" || p.find('/') != std::string::npos ? 1.0 : 1.0 / (p[0] - '0');
}

inline std::string emit(Tree& t, int node, Function& order, bool byte_order) {
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
      if (ranks[a] != ranks[b]) return ranks[a] > ranks[b];
      return !t.symmetry.empty() && t.symmetry[kids[a]] > t.symmetry[kids[b]];
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
    return result + residue_label(t, node);
  }
  t.vertices.push_back(node);
  return residue_label(t, node);
}

inline void cache_tree(Tree& t, Function& order, bool byte_order) {
  t.depth.assign(t.mono.size(), -1);
  t.sig.resize(t.mono.size());
  std::function<void(int)> visit = [&](int n) {
    if (t.depth[n] >= 0) return;
    int occupied = 0;
    std::vector<std::string> tokens;
    t.depth[n] = 0;
    for (int k : t.children[n]) {
      auto pos = t.link[k].substr(3);
      if (pos.size() == 1 && pos != "?") {
        int bit = 1 << (pos[0]-'0');
        if (occupied & bit) throw std::runtime_error("duplicated linkage positions");
        occupied |= bit;
      }
      visit(k);
      t.depth[n] = std::max(t.depth[n], t.depth[k]+1);
      tokens.push_back(t.link[k] + "->" + t.sig[k]);
    }
    t.sig[n] = residue_label(t, n);
    if (!tokens.empty()) {
      auto ids = lexical_order(tokens, order, byte_order);
      t.sig[n] += "{";
      for (int j = 0; j < (int)ids.size(); ++j) t.sig[n] += (j ? "," : "") + tokens[ids[j]];
      t.sig[n] += "}";
    }
  };
  for (int n = 0; n < (int)t.mono.size(); ++n) visit(n);
}

inline List serialize_tree(Tree& t) {
  int n = t.mono.size();
  CharacterVector monos(n), subs(n), links(t.edge_nodes.size());
  IntegerVector parents(n), edges(2*t.edge_nodes.size());
  std::vector<int> inverse(n);
  for (int j = 0; j < n; ++j) { inverse[t.vertices[j]] = j+1; monos[j] = t.mono[t.vertices[j]]; subs[j] = t.sub[t.vertices[j]]; }
  for (int j = 0; j < n; ++j) {
    int p = t.parent[t.vertices[j]]; parents[j] = p < 0 ? 0 : inverse[p];
  }
  for (int j = 0; j < (int)t.edge_nodes.size(); ++j) {
    int k = t.edge_nodes[j]; edges[2*j] = inverse[t.parent[k]]; edges[2*j+1] = inverse[k]; links[j] = t.link[k];
  }
  return List::create(_["status"]="ok", _["iupac"]=t.sequence, _["mono"]=monos,
                      _["sub"]=subs, _["parent"]=parents, _["edges"]=edges, _["linkage"]=links,
                      _["anomer"]=t.anomer, _["alditol"]=t.alditol);
}


struct Part { int root; std::vector<int> nodes, parents; std::string linkage; };
struct FloatingSub { std::string token; std::vector<int> parents; };
struct Forest { Tree tree; int root=-1; std::vector<Part> parts; std::vector<FloatingSub> subs; std::vector<int> edge_order; };

inline std::vector<int> candidates(const Part& p, int n) {
  if (!p.parents.empty()) return p.parents;
  std::vector<int> out;
  std::set<int> own(p.nodes.begin(),p.nodes.end());
  for(int i=0;i<n;++i) if(!own.count(i)) out.push_back(i);
  return out;
}
inline std::vector<int> candidates(const FloatingSub& s, int n) {
  if (!s.parents.empty()) return s.parents;
  std::vector<int> out(n); std::iota(out.begin(),out.end(),0); return out;
}
inline std::vector<int> slots(int parent, const std::vector<int>& pos) {
  std::vector<int> out; for(int p:pos) out.push_back(parent*10+p); return out;
}
inline std::vector<int> membership(const Forest& f) {
  std::vector<int> out(f.tree.mono.size(),-1);
  for(int i=0;i<(int)f.parts.size();++i) for(int n:f.parts[i].nodes) out[n]=i;
  return out;
}

inline void validate_forest(const Forest& f) {
  const Tree& t=f.tree; int n=t.mono.size(), np=f.parts.size();
  std::vector<std::vector<int>> fixed;
  for(int k:f.edge_order) {auto p=positions(t.link[k].substr(3)); if(!p.empty()) fixed.push_back(slots(t.parent[k],p));}
  for(int i=0;i<n;++i) for(auto token:split(t.sub[i],',')) {auto p=positions(token); if(!p.empty()) fixed.push_back(slots(i,p));}
  if(!matching(fixed)) throw std::runtime_error("conflicting fixed carbon positions");
  std::set<int> occupied;
  for(const auto& d:fixed) if(d.size()==1) occupied.insert(d[0]);
  std::vector<std::vector<int>> parent_domains, position_domains;
  auto add_domain = [&](const std::vector<int>& parents, const std::vector<int>& pos, bool explicit_parents) {
    if(parents.empty()) throw std::runtime_error("empty parent domain");
    std::vector<int> feasible;
    for(int parent:parents) {
      auto d=slots(parent,pos); bool free=pos.empty();
      for(int slot:d) if(!occupied.count(slot)) free=true;
      if(explicit_parents && !free) throw std::runtime_error("impossible explicit parent metadata");
      if(free) feasible.push_back(parent);
    }
    parent_domains.push_back(feasible); position_domains.push_back(pos);
  };
  for(const auto& p:f.parts) add_domain(candidates(p,n),positions(p.linkage.substr(3)),!p.parents.empty());
  for(const auto& s:f.subs) add_domain(candidates(s,n),positions(s.token),!s.parents.empty());
  auto mem=membership(f);
  // Reachability rejects disconnected candidate cycles without enumerating them.
  std::vector<bool> reaches(np,false); bool changed=true;
  while(changed) {changed=false; for(int i=0;i<np;++i) if(!reaches[i]) for(int p:parent_domains[i])
    if(mem[p]<0 || reaches[mem[p]]) {reaches[i]=true;changed=true;break;}}
  if(std::find(reaches.begin(),reaches.end(),false)!=reaches.end()) throw std::runtime_error("floating components cannot reach main");
  std::vector<int> selected(np,-2);
  std::vector<std::vector<int>> domains=fixed;
  size_t visits=0;
  std::function<bool(int)> search = [&](int i) {
    if(++visits % 4096==0) checkUserInterrupt();
    if(i==(int)parent_domains.size()) return matching(domains);
    for(int parent:parent_domains[i]) {
      bool cyclic=false;
      if(i<np) {
        selected[i]=mem[parent]; int cur=i; std::set<int> seen;
        while(cur>=0 && selected[cur]!=-2) {
          if(!seen.insert(cur).second) {cyclic=true;break;} cur=selected[cur];
        }
      }
      if(!cyclic) {
        bool constrained=!position_domains[i].empty();
        if(constrained) domains.push_back(slots(parent,position_domains[i]));
        bool ok=(!constrained || matching(domains)) && search(i+1);
        if(constrained) domains.pop_back();
        if(ok) return true;
      }
    }
    if(i<np) selected[i]=-2;
    return false;
  };
  if(!search(0)) throw std::runtime_error("no conflict-free acyclic floating assignment");
}

inline void resolve_singletons(Forest& f) {
  int n=f.tree.mono.size();
  while(true) {
    bool changed=false;
    std::vector<FloatingSub> remaining_subs;
    for(const auto& s:f.subs) {
      auto c=candidates(s,n);
      if(c.size()==1) {
        auto tokens=split(f.tree.sub[c[0]],','); tokens.push_back(s.token);
        f.tree.sub[c[0]]=collapse_subs(tokens); changed=true;
      } else remaining_subs.push_back(s);
    }
    f.subs=remaining_subs;
    std::vector<Part> remaining;
    for(const auto& p:f.parts) {
      auto c=candidates(p,n);
      if(c.size()==1) {
        f.tree.parent[p.root]=c[0]; f.tree.link[p.root]=p.linkage;
        f.tree.children[c[0]].push_back(p.root); f.edge_order.push_back(p.root); changed=true;
      } else remaining.push_back(p);
    }
    if(!changed) break;
    f.parts=remaining;
    // Only prune after localization: invalid raw explicit domains still fail
    // validation before normalization.
    std::set<int> occupied;
    for(int k:f.edge_order) {
      auto pos=positions(f.tree.link[k].substr(3));
      if(pos.size()==1) occupied.insert(f.tree.parent[k]*10+pos[0]);
    }
    for(int i=0;i<n;++i) for(const auto& token:split(f.tree.sub[i],',')) {
      auto pos=positions(token);
      if(pos.size()==1) occupied.insert(i*10+pos[0]);
    }
    auto prune = [&](std::vector<int>& parents, const std::vector<int>& pos,
                     const std::vector<int>& candidates) {
      if(pos.empty()) return;
      std::vector<int> keep;
      for(int parent:candidates) {
        for(int slot:slots(parent,pos)) if(!occupied.count(slot)) {
          keep.push_back(parent); break;
        }
      }
      if(keep.empty()) throw std::runtime_error("singleton attachment empties explicit domain");
      if(keep.size()!=candidates.size()) parents=keep;
    };

    for(auto& p:f.parts) {
      p.nodes.clear(); std::vector<int> pending={p.root};
      while(!pending.empty()) {int v=pending.back();pending.pop_back();p.nodes.push_back(v);
        for(int k:f.tree.children[v]) pending.push_back(k);}
      std::sort(p.nodes.begin(),p.nodes.end());
      if(!p.parents.empty()) {
        std::vector<int> keep;
        for(int parent:p.parents) if(!std::binary_search(p.nodes.begin(),p.nodes.end(),parent)) keep.push_back(parent);
        if(keep.empty()) throw std::runtime_error("singleton attachment empties explicit domain");
        p.parents=keep;
      }
      prune(p.parents,positions(p.linkage.substr(3)),candidates(p,n));
    }
    for(auto& s:f.subs) prune(s.parents,positions(s.token),candidates(s,n));
  }
}

inline std::vector<int> symmetry_labels(Forest& f, Function& bliss) {
  const auto& t=f.tree; int n=t.mono.size(); auto mem=membership(f);
  std::vector<std::string> keys; std::vector<int> edges;
  auto connect=[&](int a,int b){edges.push_back(a+1);edges.push_back(b+1);};
  for(int i=0;i<n;++i) {
    std::string role=mem[i]<0 ? "main-node" : "floating-node";
    if(i==f.root) role="main-root";
    else if(t.parent[i]<0) role="floating-root";
    keys.push_back("residue\r"+role+"\r"+t.mono[i]+"\r"+t.sub[i]);
  }
  for(int k:f.edge_order) {int id=keys.size();keys.push_back("glycan-edge\r"+t.link[k]);connect(t.parent[k],id);connect(id,k);}
  int part_offset=keys.size();
  for(const auto& p:f.parts) keys.push_back("floating-part\r"+p.linkage);
  int sub_offset=keys.size();
  for(const auto& s:f.subs) keys.push_back("floating-substituent\r"+s.token);
  auto relation=[&](int constraint,int target,const std::string& role) {
    int id=keys.size();keys.push_back(role);connect(constraint,id);connect(id,target);
  };
  for(int i=0;i<(int)f.parts.size();++i) {
    relation(part_offset+i,f.parts[i].root,"floating-part-root");
    for(int p:candidates(f.parts[i],n)) relation(part_offset+i,p,"floating-part-parent");
  }
  for(int i=0;i<(int)f.subs.size();++i)
    for(int p:candidates(f.subs[i],n)) relation(sub_offset+i,p,"floating-substituent-parent");
  auto unique=keys;std::sort(unique.begin(),unique.end());unique.erase(std::unique(unique.begin(),unique.end()),unique.end());
  std::vector<int> colors;for(const auto& key:keys) colors.push_back(std::lower_bound(unique.begin(),unique.end(),key)-unique.begin()+1);
  IntegerVector labels=bliss(wrap(edges),wrap(colors));
  return as<std::vector<int>>(labels);
}

inline List finish_forest(Forest& f, Function& order, bool byte_order, Function& bliss) {
  Tree& t=f.tree;
  resolve_singletons(f);
  cache_tree(t,order,byte_order);
  t.vertices.clear();t.edge_nodes.clear();
  if(f.parts.empty() && f.subs.empty()) {
    t.sequence=emit(t,f.root,order,byte_order)+(t.alditol?"-ol":"")+"("+t.anomer+"-";
    return serialize_tree(t);
  }
  auto labels=symmetry_labels(f,bliss);
  int n=t.mono.size(), part_offset=n+f.edge_order.size(), sub_offset=part_offset+f.parts.size();
  t.symmetry.assign(labels.begin(),labels.begin()+n);
  auto key_candidates=[&](const std::vector<int>& candidates) {
    std::vector<int> lab;for(int v:candidates) lab.push_back(labels[v]);std::sort(lab.begin(),lab.end());return int_join(lab);
  };
  std::vector<std::string> keys;
  std::vector<std::vector<int>> part_vertices,part_edges;
  std::vector<std::string> sequences;
  for(const auto& p:f.parts) {
    t.vertices.clear();t.edge_nodes.clear();
    std::string seq=emit(t,p.root,order,byte_order);
    part_vertices.push_back(t.vertices);part_edges.push_back(t.edge_nodes);sequences.push_back(seq);
    keys.push_back(seq+"\r"+p.linkage+"\r"+(p.parents.empty()?"all":"explicit")+"\r"+key_candidates(candidates(p,n)));
  }
  auto sorted_constraints=[&](const std::vector<std::string>& keys,int offset) {
    auto ids=lexical_order(keys,order,byte_order);
    // R orders by key then BLISS label; exact equal keys are the relevant ties.
    size_t begin=0;
    while(begin<ids.size()) {size_t end=begin+1;while(end<ids.size() && keys[ids[begin]]==keys[ids[end]]) ++end;
      std::stable_sort(ids.begin()+begin,ids.begin()+end,[&](int a,int b){return labels[offset+a]<labels[offset+b];});begin=end;}
    return ids;
  };
  auto part_ids=sorted_constraints(keys,part_offset);
  keys.clear();for(const auto& s:f.subs) keys.push_back(s.token+"\r"+(s.parents.empty()?"all":"explicit")+"\r"+key_candidates(candidates(s,n)));
  auto sub_ids=sorted_constraints(keys,sub_offset);
  t.vertices.clear();t.edge_nodes.clear();
  for(int id:part_ids) {t.vertices.insert(t.vertices.end(),part_vertices[id].begin(),part_vertices[id].end());t.edge_nodes.insert(t.edge_nodes.end(),part_edges[id].begin(),part_edges[id].end());}
  std::string main=emit(t,f.root,order,byte_order)+(t.alditol?"-ol":"")+"("+t.anomer+"-";
  std::vector<int> inverse(n);for(int i=0;i<n;++i) inverse[t.vertices[i]]=i+1;
  auto remap=[&](const std::vector<int>& parents) {std::vector<int> out;for(int p:parents) out.push_back(inverse[p]);std::sort(out.begin(),out.end());return out;};
  List parts(part_ids.size()),subs(sub_ids.size());std::string sequence;
  for(int j=0;j<(int)sub_ids.size();++j) {
    const auto& s=f.subs[sub_ids[j]];auto parents=remap(s.parents);
    subs[j]=List::create(_["substituent"]=s.token,_["parents"]=wrap(parents));
    sequence+="{"+s.token+(parents.empty()?"":"|"+int_join(parents))+"}";
  }
  for(int j=0;j<(int)part_ids.size();++j) {
    int id=part_ids[j];const auto& p=f.parts[id];auto parents=remap(p.parents);
    std::vector<int> nodes;for(int v:part_vertices[id]) nodes.push_back(inverse[v]);
    parts[j]=List::create(_["root"]=inverse[p.root],_["nodes"]=wrap(nodes),_["linkage"]=p.linkage,_["parents"]=wrap(parents));
    sequence+="{"+sequences[id]+"("+p.linkage+")"+(parents.empty()?"":"|"+int_join(parents))+"}";
  }
  t.sequence=sequence+main;
  List result=serialize_tree(t);
  if(parts.size()) result["floating_parts"]=parts;
  if(subs.size()) result["floating_substituents"]=subs;
  return result;
}


} // namespace glyrepr_compact

#endif
