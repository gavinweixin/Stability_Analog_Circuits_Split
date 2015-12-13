#include "Circuit_F2_1.h"

Circuit_F2_1 :: Circuit_F2_1()
{
    p.resize(SIZE_PARM());
}

Circuit_F2_1 :: Circuit_F2_1(const Circuit_F2_1 &orig) : Circuit(orig) { }

Circuit_F2_1 :: Circuit_F2_1(const vector<Interval>& orig) : Circuit(orig) { }

Circuit_F2_1 :: ~Circuit_F2_1() { }

Circuit_F2_1& Circuit_F2_1 :: operator = (const Circuit_F2_1& orig)
{
    if (this == &orig)
        return *this;
    p = orig.p;
    return *this;
}

size_t Circuit_F2_1 :: SIZE_PARM() const
{
    return SIZE_PARM_F2_1;
}

size_t Circuit_F2_1 :: SIZE_RT() const
{
    return SIZE_RT_F2_1;
}

vector<Interval> Circuit_F2_1 :: coef_cal() const
{
    vector<AAF> p(SIZE_PARM());
    for (size_t i=0; i<SIZE_PARM(); i++)
        p[i] = interval(this->p[i].lower(),this->p[i].upper());
    AAF coef_origin[SIZE_RT()];
    coef_origin[0] = p[1]*p[2]+p[0]*p[2];
    coef_origin[1] = p[5]*p[0]*p[1]*p[2] + p[6]*p[4]*p[1]*p[2] + p[6]*p[0]*p[1]*p[2] - p[6]*p[0]*p[4]*p[3];
    coef_origin[2] = p[6]*p[5]*p[0]*p[4]*p[1]*p[2];

    AAF coef_shift[SIZE_RT()];
    coef_shift[0] = coef_origin[2]*shiftDis*shiftDis-coef_origin[1]*shiftDis+coef_origin[0];
    coef_shift[1] = coef_origin[1]-coef_origin[2]*2.*shiftDis;
    coef_shift[2] = coef_origin[2];

    vector<Interval> coef(SIZE_RT());
    for (size_t i=0; i<SIZE_RT(); i++)
        coef[i] = Interval(coef_shift[i].convert().left(),coef_shift[i].convert().right());

    return coef;
}

vector<Interval> Circuit_F2_1 :: RouthTable() const
{
    return coef_cal();
}

vector< vector<Interval> > Circuit_F2_1 :: Jacobi_cal() const
{
    vector< vector<Interval> > temp;
    temp.resize(SIZE_RT());
    for (size_t i=0; i<SIZE_RT(); i++)
        temp[i].resize(SIZE_PARM());
    double d = shiftDis;
    temp[0][0] = p[5]*p[6]*p[1]*p[2]*p[4]*d*d+(p[5]*p[1]*p[2]+p[6]*p[1]*p[2]-p[6]*p[3]*p[4])*d+p[2];
    temp[0][1] = p[5]*p[6]*p[0]*p[2]*p[4]*d*d+(p[5]*p[0]*p[2]+p[6]*p[0]*p[2]+p[6]*p[2]*p[4])*d+p[2];
    temp[0][2] = p[5]*p[6]*p[0]*p[1]*p[4]*d*d+(p[5]*p[0]*p[1]+p[6]*p[0]*p[1]+p[6]*p[1]*p[4])*d+p[0]+p[1];
    temp[0][3] = 0.-p[6]*p[0]*p[4]*d;
    temp[0][4] = p[5]*p[6]*p[0]*p[1]*p[2]*d*d+(p[6]*p[1]*p[2]-p[6]*p[0]*p[3])*d;
    temp[0][5] = p[6]*p[0]*p[1]*p[2]*p[4]*d*d+p[0]*p[1]*p[2]*d;
    temp[0][6] = p[5]*p[0]*p[1]*p[2]*p[4]*d*d+(p[0]*p[1]*p[2]-p[0]*p[3]*p[4]+p[1]*p[2]*p[4])*d;

    temp[1][0] = p[5]*p[1]*p[2]+p[6]*p[1]*p[2]-p[6]*p[3]*p[4]+2.*p[5]*p[6]*p[1]*p[2]*p[4]*d;
    temp[1][1] = p[5]*p[0]*p[2]+p[6]*p[0]*p[2]+p[6]*p[2]*p[4]+2.*p[5]*p[6]*p[0]*p[2]*p[4]*d;
    temp[1][2] = p[5]*p[0]*p[1]+p[6]*p[0]*p[1]+p[6]*p[1]*p[4]+2.*p[5]*p[6]*p[0]*p[1]*p[4]*d;
    temp[1][3] = 0.-p[6]*p[0]*p[4];
    temp[1][4] = p[6]*p[1]*p[2]-p[6]*p[0]*p[3]+2.*p[5]*p[6]*p[0]*p[1]*p[2]*d;
    temp[1][5] = p[0]*p[1]*p[2]+2.*p[6]*p[0]*p[1]*p[2]*p[4]*d;
    temp[1][6] = p[0]*p[1]*p[2]-p[0]*p[3]*p[4]+p[1]*p[2]*p[4]+2.*p[5]*p[0]*p[1]*p[2]*p[4]*d;

    temp[2][0] = p[5]*p[6]*p[1]*p[2]*p[4];
    temp[2][1] = p[5]*p[6]*p[0]*p[2]*p[4];
    temp[2][2] = p[5]*p[6]*p[0]*p[1]*p[4];
    temp[2][3] = 0.;
    temp[2][4] = p[5]*p[6]*p[0]*p[1]*p[2];
    temp[2][5] = p[6]*p[0]*p[1]*p[2]*p[4];
    temp[2][6] = p[5]*p[0]*p[1]*p[2]*p[4];

    return temp;
}
