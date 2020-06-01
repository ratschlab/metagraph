#include "gtest_patch.hpp"

// Google Test doesn't define typeinfo for sdsl::uint128_t, so we do it for them
namespace testing {
namespace internal {

template <>
std::string GetTypeName<unsigned __int128>() { return "uint128_t"; }

} // namespace internal
} // namespace testing
