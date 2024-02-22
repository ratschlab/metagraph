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