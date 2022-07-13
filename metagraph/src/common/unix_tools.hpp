#ifndef __UNIX_TOOLS_HPP__
#define __UNIX_TOOLS_HPP__

/**
 * This header file collects some useful functions to find out about
 * system resources and stuff
 */

#include <chrono>


bool stderr_to_terminal();

/**
 * Returns the current resident set size (physical memory use) measured
 * in bytes, or zero if the value cannot be determined on this OS.
 */
size_t get_curr_RSS();
/**
 * Returns the peak (maximum so far) resident set size (physical
 * memory use) measured in bytes, or zero if the value cannot be
 * determined on this OS.
 */
size_t get_peak_RSS();


size_t get_max_files_open();

int get_num_fds();

/**
 * Estimate how many threads can be used without exceeding the max number of
 * open files, if each thread needs to open `fd_per_thread` files.
 * If the limit is too low and the limit is below `fd_per_thread`, print error and exit.
 */
int check_fd_and_adjust_threads(size_t num_threads, size_t fd_per_thread);


class Timer {
  public:
    Timer() { reset(); }

    void reset() {
        time_point_ = clock_::now();
    }

    double elapsed() const {
        return std::chrono::duration<double>(clock_::now() - time_point_).count();
    }

  private:
    typedef std::chrono::high_resolution_clock clock_;

    std::chrono::time_point<clock_> time_point_;
};


#endif // __UNIX_TOOLS_HPP__
