#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include "Circuit_F2_1.h"

using namespace std;

const double epsilon = 1e-4;
double V;
ofstream fout("axis");

pair<Circuit_F2_1, Circuit_F2_1> Bisect_j (Circuit_F2_1& orig, int j, double pos)
{
    Circuit_F2_1 temp1(orig), temp2(orig);
    double mid_j;
    if (pos == -1)
        mid_j = median(orig.get_pi(j));
    else
        mid_j = pos;
    temp1.set_pi(j, Interval(temp1.get_pi(j).lower(), mid_j));
    temp2.set_pi(j, Interval(mid_j, temp2.get_pi(j).upper()));
    pair<Circuit_F2_1, Circuit_F2_1> temp;
    temp = make_pair(temp1, temp2);
    return temp;
}

pair<Circuit_F2_1, Circuit_F2_1> Jacobi (Circuit_F2_1& orig)
{
    double max_df = 0;
    int max_j = 0;
    double d;
    vector<Interval> RT = orig.RouthTable();
    vector< vector<Interval> > df2;
    Circuit_F2_1 b_subt_bc = orig.b_sub_bc();
    df2 = orig.Jacobi_cal();
    for (int i = 0; i != SIZE_RT_F2_1; ++i)
    {
        if (RT[i] > 0) continue;
        for (int j = 0; j != SIZE_PARM_F2_1; ++j)
        {
            df2[i][j] = df2[i][j] * b_subt_bc.get_pi(j);
            d = max(fabs(df2[i][j].lower()), fabs(df2[i][j].upper()));
            if (d > max_df)
            {
                max_df = d;
                max_j = j;
            }
        }
    }
    return Bisect_j(orig,max_j,-1);
}

vector<double> findZeroF2_1I(const Circuit_F2_1 &ic);

pair<Circuit_F2_1, Circuit_F2_1> CFBM (Circuit_F2_1& orig)
{
    //RouthTable here has only one entry
    double tmp, min_pos, min_v=2*SIZE_RT_F2_1;
    int min_i=0;
    bool filter[SIZE_RT_F2_1] = {};
    double wid_RT[SIZE_RT_F2_1];
    bool zeroFound = false;

    vector<double> pos(SIZE_PARM_F2_1, -1);
    #ifdef findZeroAdded
    pos = findZeroF2_1I(orig);
    #endif
    for (int i=0; i<SIZE_PARM_F2_1; i++)
        if (pos[i]!=-1)
        {
            zeroFound = true;
        }

    vector<Interval> RT = orig.RouthTable();
    {
        using namespace boost::numeric::interval_lib::compare::certain;
        for (int i=0; i<SIZE_RT_F2_1; i++)
        {
            wid_RT[i] = width(RT[i]);
            if (RT[i] > 0.) filter[i] = true;
        }
    }
    pair<Circuit_F2_1, Circuit_F2_1> children;
    vector<Interval> RT_left, RT_right;
    for (int i = 0; i != SIZE_PARM_F2_1; ++ i)
    {
        if (zeroFound && pos[i]==-1) continue;
        children = Bisect_j(orig, i, pos[i]);
        RT_left = children.first.RouthTable();
        RT_right = children.second.RouthTable();
        tmp = 0;
        for (int j=0; j<SIZE_RT_F2_1; j++)
            if (!filter[j])
            {
                tmp += width(RT_left[j])  / wid_RT[j];
                tmp += width(RT_right[j]) / wid_RT[j];
            }
        if (tmp < min_v)
        {
            min_v = tmp;
            min_i = i;
            min_pos = pos[i];
        }
    }
    fout << min_i+1 << endl;
    return Bisect_j(orig,min_i,min_pos);
}

void Judge (Circuit_F2_1& parent, vector<Circuit_F2_1>& s, vector<Circuit_F2_1>& us, vector<Circuit_F2_1>& uc)
{
    pair<Circuit_F2_1, Circuit_F2_1> children;
    int state = parent.judge();
    if (state == 1)
    {
        s.push_back(parent);
        return;
    }
    else if (state == -1)
    {
        us.push_back(parent);
        return;
    }
    else if (parent.volume_cal() / V < epsilon)
    {
        uc.push_back(parent);
        return;
    }
    else
    {
        #ifdef SplitMethod_CFBM
        children = CFBM(parent);
        #elif defined SplitMethod_Jacobi
        children = Jacobi(parent);
        #endif
        Judge(children.first,s,us,uc);
        Judge(children.second,s,us,uc);
    }
}

Circuit_F2_1 init()
{
    vector<Interval> temp_v;
    Circuit_F2_1 temp;
    temp_v.push_back(Interval(10.00e+3*0.85, 10.00e+3*1.15));
    temp_v.push_back(Interval(10.00e+3*0.85, 10.00e+3*1.15));
    temp_v.push_back(Interval(10.00e+3*0.85, 10.00e+3*1.15));
    temp_v.push_back(Interval(588.1e+3*0.85, 588.1e+3*1.15));
    temp_v.push_back(Interval(245.4   *0.85, 245.4   *1.15));
    temp_v.push_back(Interval(90.87e-9*0.85, 90.87e-9*1.15));
    temp_v.push_back(Interval(90.07e-9*0.85, 90.07e-9*1.15));

    temp = Circuit_F2_1(temp_v);
    return temp;
}

void cubePrint(vector<Circuit_F2_1> &s, vector<Circuit_F2_1> &us, vector<Circuit_F2_1> &uc)
{
    cout << "stable:" << endl;
    for (vector<Circuit_F2_1>::iterator ivec = s.begin(); ivec != s.end(); ++ ivec)
    {
        for (int i = 0; i != 7; ++ i)
            cout << "[" << (*ivec).get_pi(i).lower() << "," << (*ivec).get_pi(i).upper() << "]" << "\t";
        cout << endl;
    }
    cout << "unstable:" << endl;
    for (vector<Circuit_F2_1>::iterator ivec = us.begin(); ivec != us.end(); ++ ivec)
    {
        for (int i = 0; i != 7; ++ i)
           cout << "[" << (*ivec).get_pi(i).lower() << "," << (*ivec).get_pi(i).upper() << "]" << "\t";
        cout << endl;
    }
    cout << "uncertain:" << endl;
    for (vector<Circuit_F2_1>::iterator ivec = uc.begin(); ivec != uc.end(); ++ ivec)
    {
        for (int i = 0; i != 7; ++ i)
            cout << "[" << (*ivec).get_pi(i).lower() << "," << (*ivec).get_pi(i).upper() << "]" << "\t";
        cout << endl;
    }
}

int main()
{
    #ifdef SplitMethod_CFBM
    cout << "Split Method:\t" << "CFBM" << endl;
    #elif defined SplitMethod_Jacobi
    cout << "Split Method:\t" << "Jacobi" << endl;
    #endif

    #ifdef IC_F2_1
    cout << "Circuit:\t\t" << "F2_1" << endl;
    #endif

    vector<Circuit_F2_1> stable;
    vector<Circuit_F2_1> unstable;
    vector<Circuit_F2_1> uncertain;
    Circuit_F2_1 p = init();

    V = p.volume_cal();
    Judge(p,stable,unstable,uncertain);

    cout << "#Cube\t" << stable.size() << "\t" << unstable.size() << "\t" << uncertain.size() << endl;

    double vol_stable=0, vol_unstable=0, vol_uncertain=0;
    for (size_t i=0; i<stable.size();    i++) vol_stable    += stable[i].volume_cal();
    for (size_t i=0; i<unstable.size();  i++) vol_unstable  += unstable[i].volume_cal();
    for (size_t i=0; i<uncertain.size(); i++) vol_uncertain += uncertain[i].volume_cal();

    cout << "Vol%\t";
    cout << setprecision(6);
    cout << vol_stable/V*100 << "%\t" << vol_unstable/V*100 << "%\t" << vol_uncertain/V*100 << "%" << endl;
    fout.close();
    return 0;
}
