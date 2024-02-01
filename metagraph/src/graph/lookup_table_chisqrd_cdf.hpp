#ifndef CHI2PLOOKUP_H
#define CHI2PLOOKUP_H

namespace mtg{
// https://www.codeproject.com/Articles/432194/How-to-Calculate-the-Chi-Squared-P-Value
// table containing pvalues for chi squared distribution with degrees of freedom 1
struct Chi2PLookup
{
    static const double * pValues[];
    static const int cutoff[];
    static const int divisor;

    inline double getPValue(double statistic, int df) {
        return((statistic >= cutoff[df-1]) ? 0.0 : pValues[df-1][int(divisor * statistic)]);
    }

};
} // namespace mtg
#endif // CHI2PLOOKUP_H

