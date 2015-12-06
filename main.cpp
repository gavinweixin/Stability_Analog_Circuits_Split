#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include "Circuit_F2_1.h"

using namespace std;

static const double epsilon = 1e-4;
double V;
ofstream fout("axis");

vector<Circuit_F2_1> Bisect_j (Circuit_F2_1& orig, int j)
{
    Circuit_F2_1 temp1(orig), temp2(orig);
    double mid_j;
    mid_j = median(orig.get_pi(j));

    temp1.set_pi(j, Interval(temp1.get_pi(j).lower(), mid_j));
    temp2.set_pi(j, Interval(mid_j, temp2.get_pi(j).upper()));
    vector<Circuit_F2_1> temp;
    temp.push_back(temp1);
    temp.push_back(temp2);
    return temp;
}

//pair<Circuit_F2_1, Circuit_F2_1> Jacobi (Circuit_F2_1& orig)
//{
//    double max_df = 0;
//    int max_j = 0;
//    double d;
//    vector<Interval> RT = orig.RouthTable();
//    vector< vector<Interval> > df2;
//    Circuit_F2_1 b_subt_bc = orig.b_sub_bc();
//    df2 = orig.Jacobi_cal();
//    for (int i = 0; i != SIZE_RT_F2_1; ++i)
//    {
//        if (RT[i] > 0) continue;
//        for (int j = 0; j != SIZE_PARM_F2_1; ++j)
//        {
//            df2[i][j] = df2[i][j] * b_subt_bc.get_pi(j);
//            d = max(fabs(df2[i][j].lower()), fabs(df2[i][j].upper()));
//            if (d > max_df)
//            {
//                max_df = d;
//                max_j = j;
//            }
//        }
//    }
//    return Bisect_j(orig,max_j,-1);
//}

vector<Circuit_F2_1> findZeroF2_1All(const Circuit_F2_1 &ic);

vector<Circuit_F2_1> CFBM (Circuit_F2_1& orig)
{
    //RouthTable here has only one entry
    double tmp, min_v=2*SIZE_RT_F2_1;
    int min_i=0;
    bool filter[SIZE_RT_F2_1] = {};
    double wid_RT[SIZE_RT_F2_1];

    #ifdef findZeroAdded
        vector<Circuit_F2_1> childrenByFindZero;
        childrenByFindZero = findZeroF2_1All(orig);
        if (childrenByFindZero.size() != 0)
            return childrenByFindZero;
    #endif

    vector<Interval> RT = orig.RouthTable();
    {
        using namespace boost::numeric::interval_lib::compare::certain;
        for (size_t i=0; i<SIZE_RT_F2_1; i++)
        {
            wid_RT[i] = width(RT[i]);
            if (RT[i] > 0.) filter[i] = true;
        }
    }
    vector<Circuit_F2_1> children(2);
    vector<Interval> RT_left, RT_right;
    for (size_t i = 0; i != SIZE_PARM_F2_1; ++ i)
    {
        children = Bisect_j(orig, i);
        RT_left = children[0].RouthTable();
        RT_right = children[1].RouthTable();
        tmp = 0;
        for (size_t j=0; j<SIZE_RT_F2_1; j++)
            if (!filter[j])
            {
                tmp += width(RT_left[j])  / wid_RT[j];
                tmp += width(RT_right[j]) / wid_RT[j];
            }
        if (tmp < min_v)
        {
            min_v = tmp;
            min_i = i;
        }
    }
    fout << min_i+1 << endl;
    return Bisect_j(orig, min_i);
}

void Judge (Circuit_F2_1& parent, vector<Circuit_F2_1>& s, vector<Circuit_F2_1>& us, vector<Circuit_F2_1>& uc)
{
    vector<Circuit_F2_1> children;
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
//        Judge(children.first,s,us,uc);
//        Judge(children.second,s,us,uc);
        size_t len = children.size();
        for (size_t i=0; i<len; i++)
            Judge(children[i],s,us,uc);
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
        for (size_t i = 0; i != SIZE_PARM_F2_1; ++ i)
            cout << "[" << (*ivec).get_pi(i).lower() << "," << (*ivec).get_pi(i).upper() << "]" << "\t";
        cout << endl;
    }
    cout << "unstable:" << endl;
    for (vector<Circuit_F2_1>::iterator ivec = us.begin(); ivec != us.end(); ++ ivec)
    {
        for (size_t i = 0; i != SIZE_PARM_F2_1; ++ i)
           cout << "[" << (*ivec).get_pi(i).lower() << "," << (*ivec).get_pi(i).upper() << "]" << "\t";
        cout << endl;
    }
    cout << "uncertain:" << endl;
    for (vector<Circuit_F2_1>::iterator ivec = uc.begin(); ivec != uc.end(); ++ ivec)
    {
        for (size_t i = 0; i != SIZE_PARM_F2_1; ++ i)
            cout << "[" << (*ivec).get_pi(i).lower() << "," << (*ivec).get_pi(i).upper() << "]" << "\t";
        cout << endl;
    }
}

int main()
{
    cout << "Split Method:\t";
    #ifdef SplitMethod_CFBM
        #ifdef findZeroAdded
        cout << "CFBM with findZero" << endl;
        #else
        cout << "CFBM" << endl;
        #endif
    #elif defined SplitMethod_Jacobi
        cout << "Jacobi" << endl;
    #endif

    cout << "Circuit:\t\t";
    #ifdef IC_F2_1
        cout << "F2_1" << endl;
    #endif

    vector<Circuit_F2_1> stable;
    vector<Circuit_F2_1> unstable;
    vector<Circuit_F2_1> uncertain;
    Circuit_F2_1 p = init();

    V = p.volume_cal();
    Judge(p,stable,unstable,uncertain);

    cout << "#Cube\t" << stable.size() << "\t" << unstable.size() << "\t" << uncertain.size() << endl;

//    vector<size_t> ratio(10, 0);
//    double vol = 0.;
//    for (vector<Circuit_F2_1>::iterator ivec = stable.begin(); ivec != stable.end(); ++ ivec)
//    {
//        ratio[(ivec->volume_cal()/V)*10000]++;
//        if (ivec->volume_cal()/V > 1e-3)
//            vol += ivec->volume_cal()/V;
//    }
//    for (size_t i=0; i<10; i++)
//        cout << i/100. << "% ~ " << (i+1)/100. << "%:\t" << ratio[i] << endl;
//    cout << "sum of vol which is larger than 0.1%: " << vol*100 << "%" << endl;

    double vol_stable=0, vol_unstable=0, vol_uncertain=0;
    for (size_t i=0; i<stable.size();    i++) vol_stable    += stable[i].volume_cal();
    for (size_t i=0; i<unstable.size();  i++) vol_unstable  += unstable[i].volume_cal();
    for (size_t i=0; i<uncertain.size(); i++) vol_uncertain += uncertain[i].volume_cal();

    cout << "Vol%\t";
    cout << setprecision(6);
    cout << vol_stable/V*100 << "%\t" << vol_unstable/V*100 << "%\t" << vol_uncertain/V*100 << "%" << endl;
    fout.close();

//    Interval x(-2,2), r(-1,1), s(-1,1);
//    Interval temp = (10.+x+r)*(10.-x+s);
//    cout << temp.lower() << "," << temp.upper() << endl;
    return 0;
}
