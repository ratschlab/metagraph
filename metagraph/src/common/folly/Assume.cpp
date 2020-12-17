// Adapted from folly

#if _USE_FOLLY

#include <folly/lang/Assume.h>

#include <iostream>

namespace folly {

namespace detail {

void assume_check(bool cond) {
  if (!cond)
    std::cerr << "compiler-hint assumption fails at runtime";
}

} // namespace detail

} // namespace folly

#endif // _USE_FOLLY
