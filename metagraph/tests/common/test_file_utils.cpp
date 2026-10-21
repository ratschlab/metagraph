#include <gtest/gtest.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>

#include <sys/wait.h>
#include <unistd.h>

#include "common/utils/file_utils.hpp"


namespace {

// Regression test for the PID guard in cleanup_tmp_dir_on_exit.
//
// create_temp_dir registers an atexit handler that removes every tracked
// temp dir. Without a PID guard, a forked child that exits normally would
// run that atexit and wipe the parent's dirs out from under it — e.g. this
// is exactly what gtest's threadsafe death-test style does. The guard makes
// the atexit a no-op in any process other than the one that registered the
// dir, so forked children can't corrupt the parent's state.
TEST(FileUtils, CreateTempDir_ForkedChildExitKeepsParentDir) {
    std::filesystem::path temp_dir = utils::create_temp_dir(
            std::filesystem::temp_directory_path(), "pid_guard_test");
    ASSERT_TRUE(std::filesystem::exists(temp_dir));

    const std::filesystem::path marker = temp_dir / "marker";
    { std::ofstream out(marker); out << "present"; }
    ASSERT_TRUE(std::filesystem::exists(marker));

    pid_t child = fork();
    ASSERT_NE(child, -1) << "fork failed";
    if (child == 0) {
        // Child process. Exit normally so atexit handlers fire. Without the
        // guard, cleanup_tmp_dir_on_exit would remove temp_dir here — which
        // the parent shares via the filesystem.
        std::exit(0);
    }

    int status = 0;
    ASSERT_EQ(waitpid(child, &status, 0), child);
    ASSERT_TRUE(WIFEXITED(status)) << "child did not exit normally";
    ASSERT_EQ(WEXITSTATUS(status), 0);

    EXPECT_TRUE(std::filesystem::exists(temp_dir))
            << "forked child's atexit wiped parent's temp dir: " << temp_dir;
    EXPECT_TRUE(std::filesystem::exists(marker))
            << "forked child's atexit wiped parent's temp dir contents: " << marker;

    // Clean up this test's dir so the process-wide atexit at shutdown has
    // nothing of ours left to do.
    utils::remove_temp_dir(temp_dir);
    EXPECT_FALSE(std::filesystem::exists(temp_dir));
}

}  // namespace
