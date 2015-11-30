#include "Interval_7D.h"
#include <iostream>
#include <cmath>

Interval_7D :: Interval_7D() { }

Interval_7D :: Interval_7D(const Interval_7D &orig) : p(orig.p) { }

Interval_7D :: Interval_7D(const vector<Interval>& orig) : p(orig) { }

Interval_7D :: ~Interval_7D() { }

Interval_7D& Interval_7D :: operator = (const Interval_7D& i)
{
    p = i.p;
    return *this;
}

double Interval_7D :: volume_cal() const
{
    return (p[0].width_cal() * p[1].width_cal() * p[2].width_cal() * p[3].width_cal() * p[4].width_cal() * p[5].width_cal() * p[6].width_cal());

}

vector<Interval> Interval_7D :: coef_cal(double d) const
{
    Interval coef_origin[SIZE_RT_F2_1];
    coef_origin[0] = p[1]*p[2]+p[0]*p[2];
    coef_origin[1] = p[5]*p[0]*p[1]*p[2] + p[6]*p[4]*p[1]*p[2] + p[6]*p[0]*p[1]*p[2] - p[6]*p[0]*p[4]*p[3];
    coef_origin[2] = p[6]*p[5]*p[0]*p[4]*p[1]*p[2];

    vector<Interval> coef;
    coef.resize(SIZE_RT_F2_1);
    coef[0] = coef_origin[2]*d*d-coef_origin[1]*d+coef_origin[0];
    coef[1] = coef_origin[1]-coef_origin[2]*2*d;
    coef[2] = coef_origin[2];

    return coef;
}

vector<Interval> Interval_7D :: RouthTable(double d) const
{
    return coef_cal(d);
}

vector< vector<Interval> > Interval_7D :: Jacobi_cal(double d) const
{
    vector< vector<Interval> > temp;
    temp.resize(SIZE_RT_F2_1);
    for (int i=0; i<SIZE_RT_F2_1; i++) temp[i].resize(SIZE_PARM_F2_1);

    temp[0][0] = p[5]*p[6]*p[1]*p[2]*p[4]*d*d+(p[5]*p[1]*p[2]+p[6]*p[1]*p[2]-p[6]*p[3]*p[4])*d+p[2];
    temp[0][1] = p[5]*p[6]*p[0]*p[2]*p[4]*d*d+(p[5]*p[0]*p[2]+p[6]*p[0]*p[2]+p[6]*p[2]*p[4])*d+p[2];
    temp[0][2] = p[5]*p[6]*p[0]*p[1]*p[4]*d*d+(p[5]*p[0]*p[1]+p[6]*p[0]*p[1]+p[6]*p[1]*p[4])*d+p[0]+p[1];
    temp[0][3] = Interval(0,0)-p[6]*p[0]*p[4]*d;
    temp[0][4] = p[5]*p[6]*p[0]*p[1]*p[2]*d*d+(p[6]*p[1]*p[2]-p[6]*p[0]*p[3])*d;
    temp[0][5] = p[6]*p[0]*p[1]*p[2]*p[4]*d*d+p[0]*p[1]*p[2]*d;
    temp[0][6] = p[5]*p[0]*p[1]*p[2]*p[4]*d*d+(p[0]*p[1]*p[2]-p[0]*p[3]*p[4]+p[1]*p[2]*p[4])*d;

    temp[1][0] = p[5]*p[1]*p[2]+p[6]*p[1]*p[2]-p[6]*p[3]*p[4]+2*p[5]*p[6]*p[1]*p[2]*p[4]*d;
    temp[1][1] = p[5]*p[0]*p[2]+p[6]*p[0]*p[2]+p[6]*p[2]*p[4]+2*p[5]*p[6]*p[0]*p[2]*p[4]*d;
    temp[1][2] = p[5]*p[0]*p[1]+p[6]*p[0]*p[1]+p[6]*p[1]*p[4]+2*p[5]*p[6]*p[0]*p[1]*p[4]*d;
    temp[1][3] = Interval(0,0)-p[6]*p[0]*p[4];
    temp[1][4] = p[6]*p[1]*p[2]-p[6]*p[0]*p[3]+2*p[5]*p[6]*p[0]*p[1]*p[2]*d;
    temp[1][5] = p[0]*p[1]*p[2]+2*p[6]*p[0]*p[1]*p[2]*p[4]*d;
    temp[1][6] = p[0]*p[1]*p[2]-p[0]*p[3]*p[4]+p[1]*p[2]*p[4]+2*p[5]*p[0]*p[1]*p[2]*p[4]*d;

    temp[2][0] = p[5]*p[6]*p[1]*p[2]*p[4];
    temp[2][1] = p[5]*p[6]*p[0]*p[2]*p[4];
    temp[2][2] = p[5]*p[6]*p[0]*p[1]*p[4];
    temp[2][3] = Interval(0,0);
    temp[2][4] = p[5]*p[6]*p[0]*p[1]*p[2];
    temp[2][5] = p[6]*p[0]*p[1]*p[2]*p[4];
    temp[2][6] = p[5]*p[0]*p[1]*p[2]*p[4];

    return temp;
}

Interval Interval_7D :: get_pi(int i) const
{
    return p[i];
}

void Interval_7D :: set_pi(int i, const Interval &value)
{
    p[i] = value;
}

Interval_7D Interval_7D :: b_sub_bc() const
{
    vector<double> inf(7);
    vector<double> sup(7);
    vector<Interval> temp_v;
    Interval_7D temp;

    for (int i = 0; i != 7; ++ i)
    {
        inf[i] = p[i].get_inf() - (p[i].get_sup()+p[i].get_inf())/2;
        sup[i] = p[i].get_sup() - (p[i].get_sup()+p[i].get_inf())/2;
        temp_v.push_back(Interval(inf[i],sup[i]));
    }

    temp = Interval_7D(temp_v);

    return temp;
}

int Interval_7D :: judge() const
{
    vector<Interval> RT = RouthTable();
    int positive=0, negtive=0, straddle=0;
    for (int i=0; i<SIZE_RT_F2_1; i++)
    {
//        if (RT[i].get_inf()>0 || abs(RT[i].get_inf())<RT[i].get_sup()/1e12) positive++;
//        else if (RT[i].get_sup()<0 || abs(RT[i].get_sup())<-RT[i].get_inf()/1e12) negtive++;
        if (RT[i].lt_0()) positive++;
        else if (RT[i].st_0()) negtive++;
        else straddle++;
    }
    if (positive==SIZE_RT_F2_1) return 1;
    else if (negtive==SIZE_RT_F2_1) return -1;
    else return 0;
}
