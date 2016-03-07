#include "Circuit_F2_1.h"

Circuit_F2_1 :: Circuit_F2_1() { }

Circuit_F2_1 :: Circuit_F2_1(const Circuit_F2_1 &orig) : p(orig.p) { }

Circuit_F2_1 :: Circuit_F2_1(const vector<Interval>& orig) : p(orig) { }

Circuit_F2_1 :: ~Circuit_F2_1() { }

Circuit_F2_1& Circuit_F2_1 :: operator = (const Circuit_F2_1& i)
{
    p = i.p;
    return *this;
}

double Circuit_F2_1 :: volume_cal() const
{
    return (width(p[0]) * width(p[1]) * width(p[2]) * width(p[3]) * width(p[4]) * width(p[5]) * width(p[6]));
}

vector<Interval> Circuit_F2_1 :: coef_cal(double d) const
{
    Interval coef_origin[SIZE_RT_F2_1];
    coef_origin[0] = p[1]*p[2]+p[0]*p[2];
    coef_origin[1] = p[5]*p[0]*p[1]*p[2] + p[6]*p[4]*p[1]*p[2] + p[6]*p[0]*p[1]*p[2] - p[6]*p[0]*p[4]*p[3];
    coef_origin[2] = p[6]*p[5]*p[0]*p[4]*p[1]*p[2];

    vector<Interval> coef;
    coef.resize(SIZE_RT_F2_1);
    coef[0] = coef_origin[2]*d*d-coef_origin[1]*d+coef_origin[0];
    coef[1] = coef_origin[1]-coef_origin[2]*2.*d;
    coef[2] = coef_origin[2];

    return coef;
}

vector<Interval> Circuit_F2_1 :: RouthTable(double d) const
{
    return coef_cal(d);
}

vector< vector<Interval> > Circuit_F2_1 :: Jacobi_cal(double d) const
{
    vector< vector<Interval> > temp;
    temp.resize(SIZE_RT_F2_1);
    for (int i=0; i<SIZE_RT_F2_1; i++)
        temp[i].resize(SIZE_PARM_F2_1);

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

Interval Circuit_F2_1 :: get_pi(int i) const
{
    return p[i];
}

void Circuit_F2_1 :: set_pi(int i, const Interval &value)
{
    p[i] = value;
}

Circuit_F2_1 Circuit_F2_1 :: b_sub_bc() const
{
    vector<Interval> temp_v;
    Circuit_F2_1 temp;

    for (int i = 0; i != SIZE_PARM_F2_1; ++ i)
    {
        temp_v.push_back(p[i]-median(p[i]));
    }

    temp = Circuit_F2_1(temp_v);

    return temp;
}

int Circuit_F2_1 :: judge(double d) const
{
    vector<Interval> RT = RouthTable(d);
    int positive=0, negtive=0, straddle=0;

    {
        using namespace boost::numeric::interval_lib::compare::certain;
        for (int i=0; i<SIZE_RT_F2_1; i++)
        {
            if (RT[i].lower()>0 || abs(RT[i].lower())<RT[i].upper()/1e12) positive++;
            else if (RT[i].upper()<0 || abs(RT[i].upper())<-RT[i].lower()/1e12) negtive++;
//            if (RT[i] > 0.) positive++;
//            else if (RT[i] < 0.) negtive++;
            else straddle++;
        }
    }
    if (positive==SIZE_RT_F2_1) return 1;
    else if (negtive==SIZE_RT_F2_1) return -1;
    else return 0;
}
