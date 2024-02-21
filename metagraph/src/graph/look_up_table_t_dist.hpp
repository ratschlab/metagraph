#include <iostream>
#include <cmath>
#include <vector>

namespace mtg{
class TDistributionTable {
private:
    std::vector<std::vector<double>> ttable;
    std::vector<double> alpha_values;

public:
    TDistributionTable();

    std::pair<double,double> getCriticalValue(double alpha, int df);
    int binary_search(double alpha);
};
}