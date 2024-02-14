#include <iostream>
#include <cmath>
#include <vector>

namespace mtg{
class TDistributionTable {
private:
    std::vector<std::vector<double>> ttable;

public:
    TDistributionTable();

    double getCriticalValue(double alpha, int df);

};
}