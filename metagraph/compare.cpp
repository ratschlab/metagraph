#include "compare.hpp"
#include "dbg_succinct_libmaus.hpp"

namespace compare {

    /**
    * Given a pointer to a graph structures G1 and G2, the function compares their elements to the
    * each other. It will perform an element wise comparison of the arrays W, last and
    * F and will only check for identity. If any element differs, the function will return 
    * false and true otherwise.
    */
    bool compare(DBG_succ* G1, DBG_succ* G2) {

        bool is_same = true;
        /*uint64_t cnt01 = 0;
        uint64_t cnt02 = 0;

        size_t top = std::min(W->n, G->W->n);
        for (size_t i = 0; i < top; ++i) {
            if (i > 0 && i % 10000 == 0) {
                std::cerr << ".";
                if (i % 100000 == 0) {
                    std::cerr << i << "/" << top << " - cnt: " << cnt01 << " " << cnt02 << std::endl;
                }
            }
        }
        std::cerr << "cnt01: " << cnt01 << std::endl;
        std::cerr << "cnt02: " << cnt02 << std::endl;
        */

        // compare size
        if (G1->W->n != G2->get_size()) {
            std::cerr << "sizes of graphs differ" << std::endl;
            std::cerr << "1: " << G1->W->n << std::endl;
            std::cerr << "2: " << G2->get_size() << std::endl;
            is_same = false;
        }
        
        // compare last
        for (size_t i = 0; i < G1->W->n; ++i) {
            if (G1->get_last(i) != G2->get_last(i)) {
                std::cerr << "last differs at position " << i << std::endl;
                std::cerr << "1: last[" << i << "] = " << G1->get_last(i)  << std::endl;
                std::cerr << "2: last[" << i << "] = " << G2->get_last(i) << std::endl;
                is_same = false;
                break;
            }
        }

        // compare W
        for (size_t i = 0; i < G1->W->n; ++i) {
            if (G1->get_W(i) != G2->get_W(i)) {
                std::cerr << "W differs at position " << i << std::endl;
                std::cerr << "1: W[" << i << "] = " << G1->get_W(i)  << std::endl;
                std::cerr << "2: W[" << i << "] = " << G2->get_W(i) << std::endl;
                is_same = false;
                break;
            }
        }

        // compare F
        for (size_t i = 0; i < G1->F.size(); ++i) {
            if (G1->get_F(i) != G2->get_F(i)) {
                std::cerr << "F differs at position " << i << std::endl;
                std::cerr << "1: F[" << i << "] = " << G1->get_F(i) << std::endl;
                std::cerr << "2: F[" << i << "] = " << G2->get_F(i) << std::endl;
                is_same = false;
                break;
            }
        }

        return is_same;

    }


}
