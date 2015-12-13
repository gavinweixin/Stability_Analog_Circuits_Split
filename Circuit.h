#ifndef CIRCUIT_H
#define CIRCUIT_H

#include <vector>
#include <aa.h>
#include <boost/numeric/interval.hpp>

typedef boost::numeric::interval<double> Interval;

using namespace std;

const double epsilon = 1e-4;
const double SHIFTD = 0;    //change the shift distance of imag axis here
#define IC_F2_1
#define SplitMethod_CFBM
//#define findZeroAdded

class Circuit
{
public:
    Circuit();
};


#endif // CIRCUIT_H
