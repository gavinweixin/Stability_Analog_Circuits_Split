#ifndef _CIRCUIT_F2_1_H_
#define _CIRCUIT_F2_1_H_

#include <Circuit.h>

class Circuit_F2_1
{
public:

    static const size_t SIZE_PARM = 7;
    static const size_t SIZE_RT = 3;
    Circuit_F2_1();
    Circuit_F2_1(const Circuit_F2_1& orig);
    Circuit_F2_1(const vector<Interval>& orig);
    ~Circuit_F2_1();
    Circuit_F2_1& operator = (const Circuit_F2_1& i);
    double volume_cal () const;
    vector<Interval> coef_cal(double=SHIFTD) const;
    vector<Interval> RouthTable(double=SHIFTD) const;
    vector< vector<Interval> > Jacobi_cal(double=SHIFTD) const;
    Interval get_pi(int i) const;
    void set_pi(int i, const Interval &value);
    Circuit_F2_1 b_sub_bc() const;
    int judge(double=SHIFTD) const;
private:
    vector<Interval> p;
};

// main algorithm
pair<Circuit_F2_1, Circuit_F2_1> Bisect_j (Circuit_F2_1& orig, int j, double pos);
pair<Circuit_F2_1, Circuit_F2_1> Jacobi (Circuit_F2_1& orig);
pair<Circuit_F2_1, Circuit_F2_1> CFBM (Circuit_F2_1& orig, double d=0);
void Judge (double,Circuit_F2_1& parent, vector<Circuit_F2_1>& s, vector<Circuit_F2_1>& us, vector<Circuit_F2_1>& uc, double d=0);

// input & output
Circuit_F2_1 init();
void normalOut(double, Circuit_F2_1 p);
void distribution(vector<Circuit_F2_1> &stable);
void shiftD(double, Circuit_F2_1 p);
void cubePrint(vector<Circuit_F2_1> &s, vector<Circuit_F2_1> &us, vector<Circuit_F2_1> &uc);

// alternate algo
vector<double> findZeroF2_1I(const Circuit_F2_1 &ic);
#endif
