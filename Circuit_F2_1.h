#ifndef _CIRCUIT_F2_1_H_
#define _CIRCUIT_F2_1_H_

#include <vector>
#include <ginac/ginac.h>
#include <boost/numeric/interval.hpp>

typedef boost::numeric::interval<double> Interval;

using namespace std;
using namespace GiNaC;

const int SIZE_PARM_F2_1 = 7;
const int SIZE_RT_F2_1 = 3;
#define IC_F2_1
#define SplitMethod_CFBM
//#define findZeroAdded

class Circuit_F2_1
{
public:
    Circuit_F2_1();
    Circuit_F2_1(const Circuit_F2_1& orig);
    Circuit_F2_1(const vector<Interval>& orig);
    ~Circuit_F2_1();
    Circuit_F2_1& operator = (const Circuit_F2_1& i);
    double volume_cal () const;
    vector<Interval> coef_cal(double) const;
    vector<Interval> RouthTable(double) const;
    vector< vector<Interval> > Jacobi_cal(double) const;
    Interval get_pi(int i) const;
    void set_pi(int i, const Interval &value);
    Circuit_F2_1 b_sub_bc() const;
    int judge(double) const;
private:
    vector<Interval> p;
};

vector<double> findZeroF2_1I(const Circuit_F2_1 &);

Interval bound(ex, const vector<symbol> &, const vector<Interval> &);
#endif
