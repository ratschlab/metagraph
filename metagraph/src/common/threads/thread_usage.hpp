#ifndef __THREAD_USAGE_HPP__
#define __THREAD_USAGE_HPP__

/**
 * The CPU time and the page faults of the CALLING thread, read from per-thread clocks: the
 * figures a request states of itself (timing.cpu_ms, timing.page_faults on /traverse, /resolve
 * and /pattern). Never the process's figures, which the concurrent requests share. A request's
 * work runs on the one thread of its handler (a multi-graph /pattern pair on one pool thread),
 * so the difference of two readings on that thread is the request's own.
 *
 * CPU: clock_gettime(CLOCK_THREAD_CPUTIME_ID) (POSIX; Linux and macOS have it), user and
 * system time of the thread, nanoseconds. Page faults: getrusage(RUSAGE_THREAD) (Linux only):
 * minor (served without I/O: a page already in memory, e.g. the page cache) and major (read
 * from the disk: an index page that was not resident). Where a figure cannot be read per
 * thread (no RUSAGE_THREAD: macOS) it is absent, never 0. Header-only.
 */

#include <cstdint>
#include <ctime>

#include <sys/resource.h>


namespace mtg {
namespace common {

struct ThreadUsage {
    bool has_cpu = false;
    uint64_t cpu_ns = 0;
    bool has_faults = false;
    uint64_t minor_faults = 0;
    uint64_t major_faults = 0;

    // the calling thread's figures now (since the thread started)
    static ThreadUsage now() {
        ThreadUsage u;
#ifdef CLOCK_THREAD_CPUTIME_ID
        timespec ts;
        if (clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts) == 0) {
            u.has_cpu = true;
            u.cpu_ns = static_cast<uint64_t>(ts.tv_sec) * 1'000'000'000
                     + static_cast<uint64_t>(ts.tv_nsec);
        }
#endif
#if defined(__linux__) && defined(RUSAGE_THREAD)
        rusage ru;
        if (getrusage(RUSAGE_THREAD, &ru) == 0) {
            u.has_faults = true;
            u.minor_faults = static_cast<uint64_t>(ru.ru_minflt);
            u.major_faults = static_cast<uint64_t>(ru.ru_majflt);
        }
#endif
        return u;
    }

    // nothing used, with the figures this platform reads per thread (the sum of no parts)
    static ThreadUsage zero() {
        const ThreadUsage u = now();
        return u.since(u);
    }

    // what the calling thread used since |start| (a reading of the same thread); a figure
    // either reading lacks is absent
    ThreadUsage since(const ThreadUsage &start) const {
        ThreadUsage d;
        d.has_cpu = has_cpu && start.has_cpu;
        if (d.has_cpu)
            d.cpu_ns = cpu_ns > start.cpu_ns ? cpu_ns - start.cpu_ns : 0;
        d.has_faults = has_faults && start.has_faults;
        if (d.has_faults) {
            d.minor_faults = minor_faults > start.minor_faults ? minor_faults - start.minor_faults : 0;
            d.major_faults = major_faults > start.major_faults ? major_faults - start.major_faults : 0;
        }
        return d;
    }

    // the figures of a part measured on another thread, added to these (a multi-graph
    // /pattern request's pairs): a figure is kept only while every part has it
    ThreadUsage& add(const ThreadUsage &other) {
        has_cpu = has_cpu && other.has_cpu;
        cpu_ns = has_cpu ? cpu_ns + other.cpu_ns : 0;
        has_faults = has_faults && other.has_faults;
        minor_faults = has_faults ? minor_faults + other.minor_faults : 0;
        major_faults = has_faults ? major_faults + other.major_faults : 0;
        return *this;
    }

    double cpu_ms() const { return static_cast<double>(cpu_ns) / 1e6; }
};

// A span of the calling thread's work: started at construction, read by elapsed() on the same
// thread
class ThreadUsageMeter {
  public:
    ThreadUsageMeter() : start_(ThreadUsage::now()) {}
    ThreadUsage elapsed() const { return ThreadUsage::now().since(start_); }
    void reset() { start_ = ThreadUsage::now(); }

  private:
    ThreadUsage start_;
};

} // namespace common
} // namespace mtg

#endif // __THREAD_USAGE_HPP__
