#ifndef GLYREPR_COMPACT_SYNTAX_H
#define GLYREPR_COMPACT_SYNTAX_H

#include <string>

namespace glyrepr_compact {

// These tokens have an ASCII grammar. Do not compile std::regex at DLL load
// time: its Windows collation machinery can hang under a non-C LC_COLLATE.
// Return the end of the first token, matching the former regex alternation
// order, or begin when no token starts here.
inline size_t substituent_end(const std::string& text, size_t begin) {
  size_t i = begin;
  if (i == text.size()) return begin;
  if (text[i] == '?') {
    ++i;
  } else {
    if (text[i] < '1' || text[i] > '9') return begin;
    ++i;
    while (i + 1 < text.size() && text[i] == '/' &&
           text[i + 1] >= '1' && text[i + 1] <= '9') i += 2;
  }
  static const char* const names[] = {
    "PPEtn", "PEtn", "NAc", "NGc", "Pyr", "Me", "Ac", "PC", "Gc", "P", "S", "N"
  };
  for (const char* name : names) {
    size_t length = std::char_traits<char>::length(name);
    if (text.compare(i, length, name) == 0) return i + length;
  }
  return begin;
}

inline bool is_substituent(const std::string& text) {
  return !text.empty() && substituent_end(text, 0) == text.size();
}

inline bool parent_indices_syntax(const std::string& text) {
  size_t i = 0;
  while (i < text.size()) {
    if (text[i] < '1' || text[i] > '9') return false;
    ++i;
    while (i < text.size() && text[i] >= '0' && text[i] <= '9') ++i;
    if (i == text.size()) return true;
    if (text[i++] != ',') return false;
  }
  return false;
}

} // namespace glyrepr_compact

#endif
