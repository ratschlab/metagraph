#ifndef __COMPARE_HPP__
#define __COMPARE_HPP__

#include <cstdint>

#include "dbg_succinct_libmaus.hpp"

namespace compare {

    /**
    * Given a pointer to a graph structures G1 and G2, the function compares their elements to the
    * each other. It will perform an element wise comparison of the arrays W, last and
    * F and will only check for identity. If any element differs, the function will return
    * false and true otherwise.
    */
    bool compare(DBG_succ* G1, DBG_succ* G2);

}

#endif
