#ifndef __SOURCE_REVISION_HPP__
#define __SOURCE_REVISION_HPP__


namespace mtg {
namespace common {

// The source a binary was built from, as CMakeListsSourceRevision.cmake read the checkout at
// build time (`metagraph --version` states them; GPT review 7, finding 6): the commit (40 hex
// digits, or "unknown" without git or a checkout), and "clean", "dirty" (a tracked change, a
// submodule off its commit, or an untracked file under src/, which CMake compiles) or
// "unknown". Defined in the generated source the build writes (gen/source_revision.cpp)
const char* source_commit();
const char* source_state();

} // namespace common
} // namespace mtg

#endif // __SOURCE_REVISION_HPP__
