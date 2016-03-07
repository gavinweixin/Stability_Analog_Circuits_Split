#include "Circuit_F2_1.h"

// ofstream fout("axis");

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
    for (size_t i = 0; i != orig.SIZE_RT; ++i)
    {
        if (RT[i] > 0) continue;
        for (size_t j = 0; j != orig.SIZE_PARM; ++j)
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

pair<Circuit_F2_1, Circuit_F2_1> CFBM (Circuit_F2_1& orig, double d)
{
    //RouthTable here has only one entry
    double tmp, min_pos, min_v=2*orig.SIZE_RT;
    int min_i=0;
    bool filter[orig.SIZE_RT] = {};
    double wid_RT[orig.SIZE_RT];
    bool zeroFound = false;

    vector<double> pos(orig.SIZE_PARM, -1);
    #ifdef findZeroAdded
    pos = findZeroF2_1I(orig);
    #endif
    for (size_t i=0; i<orig.SIZE_PARM; i++)
        if (pos[i]!=-1)
        {
            zeroFound = true;
        }

    vector<Interval> RT = orig.RouthTable(d);
    {
        using namespace boost::numeric::interval_lib::compare::certain;
        for (size_t i=0; i<orig.SIZE_RT; i++)
        {
            wid_RT[i] = width(RT[i]);
            if (RT[i] > 0.) filter[i] = true;
        }
    }
    pair<Circuit_F2_1, Circuit_F2_1> children;
    vector<Interval> RT_left, RT_right;
    for (size_t i = 0; i != orig.SIZE_PARM; ++ i)
    {
        if (zeroFound && pos[i]==-1) continue;
        children = Bisect_j(orig, i, pos[i]);
        RT_left = children.first.RouthTable(d);
        RT_right = children.second.RouthTable(d);
        tmp = 0;
        for (size_t j=0; j<orig.SIZE_RT; j++)
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
//    fout << min_i+1 << endl;
    return Bisect_j(orig,min_i,min_pos);
}

void Judge (double totVol, Circuit_F2_1& parent, vector<Circuit_F2_1>& s, vector<Circuit_F2_1>& us, vector<Circuit_F2_1>& uc, double d)
{
    pair<Circuit_F2_1, Circuit_F2_1> children;
    int state = parent.judge(d);
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
    else if (parent.volume_cal() / totVol < epsilon)
    {
        uc.push_back(parent);
        return;
    }
    else
    {
        #ifdef SplitMethod_CFBM
        children = CFBM(parent, d);
        #elif defined SplitMethod_Jacobi
        children = Jacobi(parent);
        #endif
        Judge(totVol,children.first,s,us,uc,d);
        Judge(totVol,children.second,s,us,uc,d);
    }
}

Circuit_F2_1 init()
{
    vector<Interval> temp_v;
    Circuit_F2_1 temp;
    temp_v.push_back(Interval(1e+3*0.85, 1e+3*1.15));
    temp_v.push_back(Interval(19e+3*0.85, 19e+3*1.15));
    temp_v.push_back(Interval(1600e-12*0.85, 1600e-12*1.15));
    temp_v.push_back(Interval(2.5133e+7*1.2*0.85, 2.5133e+7*1.2*1.15));
    temp_v.push_back(Interval(1e4*0.85, 1e4*1.15));

    temp = Circuit_F2_1(temp_v);
    return temp;
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

    Circuit_F2_1 p = init();
    double totVol = p.volume_cal();
    normalOut(totVol, p);

//    fout.close();
    return 0;
}
