#ifndef CIRCUIT_H
#define CIRCUIT_H

#include <vector>
#include <queue>
#include <aa.h>
#include <boost/numeric/interval.hpp>

typedef boost::numeric::interval<double> Interval;

using namespace std;

const double epsilon = 1e-4;
//const double SHIFTD = 0;    //change the shift distance of imag axis here
#define IC_F2_1
#define SplitMethod_CFBM
//#define findZeroAdded

class Circuit
{
public:
    Circuit();
    Circuit(const Circuit&);
    Circuit(const vector<Interval>&);
    ~Circuit();
    Circuit& operator = (const Circuit&);
    Interval get_pi(int i) const;
    void set_pi(int i, const Interval &value);
    virtual int judge(double=SHIFTD) const;
    virtual double volume_cal () const;

    virtual size_t SIZE_PARM() const = 0;
    virtual size_t SIZE_RT() const = 0;
    virtual vector<Interval> coef_cal() const = 0;
    virtual vector<Interval> RouthTable() const = 0;
    virtual vector< vector<Interval> > Jacobi_cal() const = 0;

    // main algorithm
    pair<Circuit, Circuit> Bisect_j (Circuit&, int, double);
    pair<Circuit, Circuit> Jacobi (Circuit&);
    pair<Circuit, Circuit> CFBM (Circuit&, double=0);
    void Judge (double=0);
protected:
    queue< vector<Interval> > p;
    double totVol;
    double shiftDis;
    vector< vector<Interval> > stable, unstable, uncertain;
};

// input & output
void normalOut(double, Circuit);
void distribution(vector<Circuit> &);
void shiftD(double, Circuit);
void cubePrint(Circuit&, vector<Circuit>&, vector<Circuit>&, vector<Circuit>&);

#endif // CIRCUIT_H
