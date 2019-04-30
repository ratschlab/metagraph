#ifndef __UNIX_TOOLS_HPP__
#define __UNIX_TOOLS_HPP__

/**
 * This header file collects some useful functions to find out about
 * system resources and stuff
 */

#include <chrono>


/**
 * Returns the current resident set size (physical memory use) measured
 * in bytes, or zero if the value cannot be determined on this OS.
 *
 * The code was copied and has been modified from:
 * http://nadeausoftware.com/articles/2012/07/c_c_tip_how_get_process_resident_set_size_physical_memory_use
 */
size_t get_curr_RSS();
/**
 * Returns the peak (maximum so far) resident set size (physical
 * memory use) measured in bytes, or zero if the value cannot be
 * determined on this OS.
 *
 * The code was copied and has been modified from:
 * http://nadeausoftware.com/articles/2012/07/c_c_tip_how_get_process_resident_set_size_physical_memory_use
 */
size_t get_peak_RSS();


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
