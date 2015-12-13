#include "Circuit.h"

Circuit :: Circuit() { }

Circuit :: Circuit(const Circuit &orig) : p(orig.p) { }

Circuit :: Circuit(const vector<Interval>& orig) : p(orig) { }

Circuit :: ~Circuit() { }

Circuit& Circuit :: operator = (const Circuit &orig)
{
    if (this == &orig)
        return *this;
    p = orig.p;
    return *this;
}

Interval Circuit :: get_pi(int i) const
{
    return p[i];
}

void Circuit :: set_pi(int i, const Interval &value)
{
    p[i] = value;
}

double Circuit :: volume_cal(size_t pos) const
{
    if (p[pos])
    return (width(p[0])*width(p[1])*width(p[2])*width(p[3])*width(p[4])*width(p[5])*width(p[6]));
}

int Circuit :: judge() const
{
    vector<Interval> RT = RouthTable();
    size_t positive=0, negtive=0, straddle=0;

    {
        using namespace boost::numeric::interval_lib::compare::certain;
        for (size_t i=0; i<SIZE_RT(); i++)
        {
//            if (RT[i].lower()>0 || abs(RT[i].lower())<RT[i].upper()/1e12) positive++;
//            else if (RT[i].upper()<0 || abs(RT[i].upper())<-RT[i].lower()/1e12) negtive++;
            if (RT[i] > 0.) positive++;
            else if (RT[i] < 0.) negtive++;
            else straddle++;
        }
    }
    if (positive==SIZE_RT()) return 1;
    else if (negtive==SIZE_RT()) return -1;
    else return 0;
}

//---------------------------------------------------------------------------------------
pair<Circuit, Circuit> Bisect_j (Circuit& orig, int j, double pos)
{
    Circuit temp1(orig), temp2(orig);
    double mid_j;
    if (pos == -1)
        mid_j = median(orig.get_pi(j));
    else
        mid_j = pos;
    temp1.set_pi(j, Interval(temp1.get_pi(j).lower(), mid_j));
    temp2.set_pi(j, Interval(mid_j, temp2.get_pi(j).upper()));
    pair<Circuit, Circuit> temp;
    temp = make_pair(temp1, temp2);
    return temp;
}

pair<Circuit, Circuit> Jacobi (Circuit& orig)
{
    double max_df = 0;
    int max_j = 0;
    double d;
    vector<Interval> RT = orig.RouthTable();
    vector< vector<Interval> > df2;
    Circuit b_subt_bc;

    for (size_t i = 0; i != orig.SIZE_PARM(); ++ i)
        b_subt_bc.set_pi(i, orig.get_pi(i)-median(orig.get_pi(i)));

    df2 = orig.Jacobi_cal();
    for (size_t i = 0; i != orig.SIZE_RT(); ++i)
    {
        if (RT[i] > 0) continue;
        for (size_t j = 0; j != orig.SIZE_PARM(); ++j)
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

pair<Circuit, Circuit> CFBM (Circuit& orig, double d)
{
    //RouthTable here has only one entry
    double tmp, min_pos, min_v=2*orig.SIZE_RT();
    int min_i=0;
    vector<bool> filter(orig.SIZE_RT(), false);
//    bool filter[orig.SIZE_RT()] = {};
    double wid_RT[orig.SIZE_RT()];
    bool zeroFound = false;

    vector<double> pos(orig.SIZE_PARM(), -1);
//    #ifdef findZeroAdded
//    pos = findZeroF2_1I(orig);
//    #endif
    for (size_t i=0; i<orig.SIZE_PARM(); i++)
        if (pos[i]!=-1)
        {
            zeroFound = true;
        }

    vector<Interval> RT = orig.RouthTable(d);
    {
        using namespace boost::numeric::interval_lib::compare::certain;
        for (size_t i=0; i<orig.SIZE_RT(); i++)
        {
            wid_RT[i] = width(RT[i]);
            if (RT[i] > 0.) filter[i] = true;
        }
    }
    pair<Circuit, Circuit> children;
    vector<Interval> RT_left, RT_right;
    for (size_t i = 0; i != orig.SIZE_PARM(); ++ i)
    {
        if (zeroFound && pos[i]==-1) continue;
        children = Bisect_j(orig, i, pos[i]);
        RT_left = children.first.RouthTable(d);
        RT_right = children.second.RouthTable(d);
        tmp = 0;
        for (size_t j=0; j<orig.SIZE_RT(); j++)
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
    return Bisect_j(orig,min_i,min_pos);
}

void Circuit :: Judge ()
{
    if (p.empty())
        return;

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
