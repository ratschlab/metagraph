#ifndef CHI2PLOOKUP_H
#define CHI2PLOOKUP_H

namespace mtg{
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

