#ifndef __TLOOKUP_HPP__
#define __TLOOKUP_HPP__

#include <iostream>
#include <cmath>
#include <vector>
#include <boost/math/distributions/students_t.hpp>

namespace mtg{
class TDistributionTable {
private:

public:
    TDistributionTable();

    double getCriticalValue(double alpha, int df);
};
}

#endif // __TLOOKUP_HPP__