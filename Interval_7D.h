#include "Interval.h"
#include <vector>
using namespace std;

#ifndef _INTERVAL_7D_H_
#define _INTERVAL_7D_H_

const int SIZE_PARM_F2_1 = 7;
const int SIZE_RT_F2_1 = 3;
const double SHIFTD = 0;    //change the shift distance of imag axis here

class Interval_7D
{
public:
    Interval_7D();
    Interval_7D(const Interval_7D& orig);
    Interval_7D(const vector<Interval>& orig);
    ~Interval_7D();
    Interval_7D& operator = (const Interval_7D& i);
    double volume_cal () const;
    Interval* RouthTable(double=SHIFTD) const;
    vector<Interval> J2_cal() const;
    Interval get_pi(int i) const;
    void set_pi(int i, const Interval &value);
    Interval_7D b_sub_bc() const;
    int judge() const;
private:
    vector<Interval> p;
};

#endif
