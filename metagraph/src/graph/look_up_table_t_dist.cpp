#include "look_up_table_t_dist.hpp"

namespace mtg{

TDistributionTable::TDistributionTable(){};
    // populate t table

double TDistributionTable::getCriticalValue(double alpha, int df){
    // Create Student's t-distribution object
    boost::math::students_t dist(df);

    // Get critical value for two-tailed test
    double critical_value = boost::math::quantile(dist, 1 - alpha);
    return critical_value;
}   

}; // namespace mtg