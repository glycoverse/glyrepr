#include "compact-core.h"

namespace glyrepr_compact {

Tree parse_tree(std::string s, const std::unordered_map<std::string, int>& known,
                const std::unordered_map<std::string, std::string>& configurations) {
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
      if (link == "?-?") link = "?" "?-?";
      // Normalize any ambiguous acceptor containing '?' to unknown, like R.
      if (link.size() >= 4 && link[2] == '-' && link.substr(3).find('?') != std::string::npos) {
        static const std::regex unknown("^[ab?][12?]-([1-9]/)*\\?(/[1-9])*$");
        if (std::regex_match(link, unknown)) link.replace(3, std::string::npos, "?");
      }
      if (!std::regex_match(link, linkage)) throw std::runtime_error("invalid linkage");
    }
    int end = cursor;
    while (cursor >= 0 && s[cursor] != '[' && s[cursor] != ']' && s[cursor] != ')') --cursor;
    std::string mono = s.substr(cursor+1, end-cursor);
    if (!std::regex_match(mono, mono_pattern)) throw std::runtime_error("invalid residue syntax");
    auto residue = extract_residue(mono, known, configurations);
    mono = residue.first;
    auto it = known.find(mono);
    if ((current == -1) != link.empty()) throw std::runtime_error("missing or misplaced linkage");
    int id = t.mono.size();
    if (id >= 2048) throw Unsupported("native size guard");
    t.mono.push_back(mono); t.sub.push_back(residue.second); t.link.push_back(link); t.parent.push_back(current);
    t.children.emplace_back(); t.depth.push_back(0); t.sig.emplace_back();
    if (current >= 0) t.children[current].push_back(id);
    else if (t.anomer.empty()) t.anomer = "?" + std::to_string(it->second);
    current = id;
    if (!branch_has_node.empty()) branch_has_node.back() = true;
  }
  if (!stack.empty() || t.mono.empty()) throw std::runtime_error("unbalanced branch or empty tree");
  return t;
}

std::vector<int> parse_parents(const std::string& text) {
  static const std::regex pattern("^[1-9][0-9]*(,[1-9][0-9]*)*$");
  if(!std::regex_match(text,pattern)) throw std::runtime_error("invalid parent indices");
  std::vector<int> out;
  for(auto token:split(text,',')) {long long p=std::stoll(token);if(p>INT_MAX) throw std::runtime_error("parent index overflow");out.push_back((int)p-1);}
  std::set<int> unique(out.begin(),out.end());if(unique.size()!=out.size()) throw std::runtime_error("duplicate parent indices");
  return out;
}

List parse_complete(std::string s,const std::unordered_map<std::string,int>& known,
    const std::unordered_map<std::string,std::string>& configurations,Function& order,bool byte_order,Function& bliss) {
  // Same two linkage normalizations as the public string entrypoint.
  size_t p=0;while((p=s.find("(?-?)",p))!=std::string::npos) {s.replace(p,5,"(?" "?-?)");p+=6;}
  static const std::regex unknown("\\(([ab?][12?])-([1-9]/)*\\?(/[1-9])*\\)");
  s=std::regex_replace(s,unknown,"($1-?)");
  if(s.empty() || s[0]!='{') {
    Tree t=parse_tree(s,known,configurations);cache_tree(t,order,byte_order);
    t.sequence=emit(t,0,order,byte_order)+(t.alditol?"-ol":"")+"("+t.anomer+"-";
    return serialize_tree(t);
  }
  Forest f;std::vector<int> source_map;
  auto append_component=[&](Tree& component) {
    cache_tree(component,order,byte_order);emit(component,0,order,byte_order);
    int n=component.mono.size(),offset=f.tree.mono.size();
    if (n + offset > 2048) throw Unsupported("native size guard");
    std::vector<int> inverse(n);
    for(int j=0;j<n;++j) inverse[component.vertices[j]]=offset+j;
    for(int id:component.vertices) {
      f.tree.mono.push_back(component.mono[id]);f.tree.sub.push_back(component.sub[id]);f.tree.link.push_back(component.link[id]);
      f.tree.parent.push_back(component.parent[id]<0?-1:inverse[component.parent[id]]);f.tree.children.emplace_back();
    }
    for(int k:component.edge_nodes) {int id=inverse[k];f.tree.children[f.tree.parent[id]].push_back(id);f.edge_order.push_back(id);}
    for(int i=n-1;i>=0;--i) source_map.push_back(inverse[i]);
    return inverse[0];
  };
  while(!s.empty() && s[0]=='{') {
    size_t end=s.find('}');if(end==std::string::npos) throw std::runtime_error("unclosed floating block");
    std::string content=s.substr(1,end-1);s.erase(0,end+1);
    if(content.find('{')!=std::string::npos) throw std::runtime_error("nested floating block");
    auto fields=split(content,'|');if(fields.empty() || fields.size()>2 || fields[0].empty()) throw std::runtime_error("invalid floating block");
    std::vector<int> parents=fields.size()==2?parse_parents(fields[1]):std::vector<int>();
    if(std::regex_match(fields[0],sub_regex)) {f.subs.push_back({normalize_sub(fields[0]),parents});continue;}
    auto sequence=fields[0];size_t start=sequence.rfind('(');
    if(start==std::string::npos || sequence.back()!=')') throw std::runtime_error("floating part requires linkage");
    std::string linkage=sequence.substr(start+1,sequence.size()-start-2);
    static const std::regex link_pattern("^[ab?][12?]-([1-9](/[1-9])*|\\?)$");
    if(!std::regex_match(linkage,link_pattern)) throw std::runtime_error("invalid floating linkage");
    Tree component=parse_tree(sequence.substr(0,start)+"("+linkage.substr(0,2)+"-",known,configurations);
    if(component.alditol) throw std::runtime_error("floating alditol");
    int offset=f.tree.mono.size(),root=append_component(component);
    std::vector<int> nodes(component.mono.size());std::iota(nodes.begin(),nodes.end(),offset);
    f.parts.push_back({root,nodes,parents,linkage});
  }
  if(s.empty() || s.find_first_of("{}")!=std::string::npos) throw std::runtime_error("missing or malformed main tree");
  Tree main=parse_tree(s,known,configurations);f.root=append_component(main);f.tree.anomer=main.anomer;f.tree.alditol=main.alditol;
  auto remap=[&](std::vector<int>& parents) {
    for(int& p:parents) {if(p>=(int)source_map.size()) throw std::runtime_error("parent outside complete structure");p=source_map[p];}
    std::sort(parents.begin(),parents.end());
  };
  for(auto& part:f.parts) {
    remap(part.parents);
    for(int p:part.parents) if(std::find(part.nodes.begin(),part.nodes.end(),p)!=part.nodes.end()) throw std::runtime_error("self parent metadata");
  }
  for(auto& sub:f.subs) remap(sub.parents);
  validate_forest(f);
  return finish_forest(f,order,byte_order,bliss);
}

} // namespace glyrepr_compact

// [[Rcpp::export(name = ".compact_parse_native")]]
List compact_parse_native(CharacterVector strings,CharacterVector residues,IntegerVector anomer_positions,
    CharacterVector configuration_names,CharacterVector configuration_values,Function order,Function bliss,bool byte_order=false) {
  std::unordered_map<std::string,int> known;
  std::unordered_map<std::string,std::string> configurations;
  for(int i=0;i<residues.size();++i) known.emplace(as<std::string>(residues[i]),anomer_positions[i]);
  for(int i=0;i<configuration_names.size();++i) configurations.emplace(as<std::string>(configuration_names[i]),as<std::string>(configuration_values[i]));
  List out(strings.size());
  for(int i=0;i<strings.size();++i) {
    if(i%128==0) checkUserInterrupt();
    if(strings[i]==NA_STRING) {out[i]=List::create(_["status"]="missing");continue;}
    try {out[i]=glyrepr_compact::parse_complete(as<std::string>(strings[i]),known,configurations,order,byte_order,bliss);}
    catch(const glyrepr_compact::Unsupported& e) {out[i]=List::create(_["status"]="unsupported",_["reason"]=e.what());}
    catch(const std::exception& e) {out[i]=List::create(_["status"]="error",_["reason"]=e.what());}
  }
  return out;
}
