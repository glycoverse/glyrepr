#ifndef GLYREPR_COMPACT_PROGRESS_H
#define GLYREPR_COMPACT_PROGRESS_H
#include <Rcpp.h>
#include <chrono>

// Invoked outside element recovery: callback failures must abort the operation.
class CompactProgress {
  SEXP callback;
  std::chrono::steady_clock::time_point last;
public:
  explicit CompactProgress(SEXP callback_) : callback(callback_) {
    if (callback != R_NilValue) last = std::chrono::steady_clock::now();
  }
  void update(int current, bool final = false) {
    if (callback == R_NilValue) return;
    auto now = std::chrono::steady_clock::now();
    if (final || std::chrono::duration<double>(now - last).count() >= 0.1) {
      Rcpp::Function report(callback);
      report(current);
      last = now;
    }
  }
};
#endif
