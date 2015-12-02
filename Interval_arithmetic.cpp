#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include "Interval_7D.h"

using namespace std;

const double epsilon = 0.0001;
double V;
ofstream fout("axis");

pair<Interval_7D, Interval_7D> Bisect_j (Interval_7D& orig, int j, double pos)
{
    Interval_7D temp1(orig), temp2(orig);
    double mid_j;
    if (pos == -1)
        mid_j = (orig.get_pi(j).get_inf() + orig.get_pi(j).get_sup()) / 2;
    else
        mid_j = pos;
    temp1.set_pi(j, Interval(temp1.get_pi(j).get_inf(), mid_j));
    temp2.set_pi(j, Interval(mid_j, temp2.get_pi(j).get_sup()));
    pair<Interval_7D, Interval_7D> temp;
    temp = make_pair(temp1, temp2);
    return temp;
}

pair<Interval_7D, Interval_7D> Jacobi (Interval_7D& orig)
{
    double max_df = 0;
    int max_j = 0;
    double d;
    vector<Interval> RT = orig.RouthTable();
    vector< vector<Interval> > df2;
    Interval_7D b_subt_bc = orig.b_sub_bc();
    df2 = orig.Jacobi_cal();
    for (int i = 0; i != SIZE_RT_F2_1; ++i)
    {
        if (RT[i].lt_0()) continue;
        for (int j = 0; j != SIZE_PARM_F2_1; ++j)
        {
            df2[i][j] = df2[i][j] * b_subt_bc.get_pi(j);
            d = max(fabs(df2[i][j].get_inf()), fabs(df2[i][j].get_sup()));
            if (d > max_df)
            {
                max_df = d;
                max_j = j;
            }
        }
    }
    return Bisect_j(orig,max_j,-1);
}

vector<double> findZeroF2_1I(const Interval_7D &ic);

pair<Interval_7D, Interval_7D> CFBM (Interval_7D& orig)
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
    for (int i=0; i<SIZE_RT_F2_1; i++)
    {
        wid_RT[i] = RT[i].width_cal();
        if (RT[i].lt_0()) filter[i] = true;
    }

    pair<Interval_7D, Interval_7D> children;
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
                tmp += RT_left[j].width_cal()  / wid_RT[j];
                tmp += RT_right[j].width_cal() / wid_RT[j];
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

void Judge (Interval_7D& parent, vector<Interval_7D>& s, vector<Interval_7D>& us, vector<Interval_7D>& uc)
{
    pair<Interval_7D, Interval_7D> children;
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

Interval_7D init()
{
    vector<Interval> temp_v;
    Interval_7D temp;
    temp_v.push_back(Interval(10.00e+3*0.85, 10.00e+3*1.15));
    temp_v.push_back(Interval(10.00e+3*0.85, 10.00e+3*1.15));
    temp_v.push_back(Interval(10.00e+3*0.85, 10.00e+3*1.15));
    temp_v.push_back(Interval(588.1e+3*0.85, 588.1e+3*1.15));
    temp_v.push_back(Interval(245.4   *0.85, 245.4   *1.15));
    temp_v.push_back(Interval(90.87e-9*0.85, 90.87e-9*1.15));
    temp_v.push_back(Interval(90.07e-9*0.85, 90.07e-9*1.15));

    temp = Interval_7D(temp_v);
    return temp;
}

void cubePrint(vector<Interval_7D> &s, vector<Interval_7D> &us, vector<Interval_7D> &uc)
{
    cout << "stable:" << endl;
    for (vector<Interval_7D>::iterator ivec = s.begin(); ivec != s.end(); ++ ivec)
    {
        for (int i = 0; i != 7; ++ i)
            cout << "[" << (*ivec).get_pi(i).get_inf() << "," << (*ivec).get_pi(i).get_sup() << "]" << "\t";
        cout << endl;
    }
    cout << "unstable:" << endl;
    for (vector<Interval_7D>::iterator ivec = us.begin(); ivec != us.end(); ++ ivec)
    {
        for (int i = 0; i != 7; ++ i)
           cout << "[" << (*ivec).get_pi(i).get_inf() << "," << (*ivec).get_pi(i).get_sup() << "]" << "\t";
        cout << endl;
    }
    cout << "uncertain:" << endl;
    for (vector<Interval_7D>::iterator ivec = uc.begin(); ivec != uc.end(); ++ ivec)
    {
        for (int i = 0; i != 7; ++ i)
            cout << "[" << (*ivec).get_pi(i).get_inf() << "," << (*ivec).get_pi(i).get_sup() << "]" << "\t";
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

    #ifdef Circuit_F2_1
    cout << "Circuit:\t\t" << "F2_1" << endl;
    #endif

    vector<Interval_7D> stable;
    vector<Interval_7D> unstable;
    vector<Interval_7D> uncertain;
    Interval_7D p = init();

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
