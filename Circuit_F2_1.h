#ifndef _CIRCUIT_F2_1_H_
#define _CIRCUIT_F2_1_H_

#include <Circuit.h>

class Circuit_F2_1 : public Circuit
{
public:
    Circuit_F2_1();
    Circuit_F2_1(const Circuit_F2_1&);
    Circuit_F2_1(const vector<Interval>&);
    ~Circuit_F2_1();
    Circuit_F2_1& operator = (const Circuit_F2_1&);

    size_t SIZE_PARM() const;
    size_t SIZE_RT() const;
    vector<Interval> coef_cal(double=SHIFTD) const;
    vector<Interval> RouthTable(double=SHIFTD) const;
    vector< vector<Interval> > Jacobi_cal(double=SHIFTD) const;

    static const size_t SIZE_PARM_F2_1 = 7;
    static const size_t SIZE_RT_F2_1 = 3;
};

Circuit_F2_1 init();
// alternate algo
vector<double> findZeroF2_1I(const Circuit_F2_1 &ic);
#endif
