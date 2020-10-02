#include <cstdint>

// work around clang-related bug in GMock
namespace testing {
namespace internal {
using UInt64 = uint64_t;
using Int64 = int64_t;
using Int32 = int32_t;
} // namespace internal
} // namespace testing
