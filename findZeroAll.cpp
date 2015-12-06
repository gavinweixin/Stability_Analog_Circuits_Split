#include <Circuit_F2_1.h>
#include <cmath>

//bounds of interval not included
bool insideInterval(const Interval &intval, const double &value)
{
    double eps = min(median(intval)/10000, 1e-8);
//    double eps = width(intval)/10;  //to be modified
    return ((value-intval.lower()>eps) && (intval.upper()-value>eps));
}

bool newtonMid(Circuit_F2_1 ic, int num, Interval range, double &zero)
{
    const size_t MAXITEA = 20;
    vector<Interval> f;
    Interval df;
    double itr;
    double eps = median(range)/10000;  //set the accuracy
    bool found = false;

    for (size_t i=0; i<SIZE_PARM_F2_1; i++)
        ic.set_pi(i, Interval(median(ic.get_pi(i))));
    ic.set_pi(num, Interval(median(range)));

    for (size_t i=0; i<MAXITEA; i++)
    {
        f = ic.RouthTable();
        df = (ic.Jacobi_cal())[1][num];    //to be modified
        itr = median(ic.get_pi(num))-median(f[1])/median(df);
        if (abs(itr-median(ic.get_pi(num)))<eps)
        {
            found = true;
            break;
        }
        ic.set_pi(num, Interval(itr));
    }
    zero = itr;
    if (!insideInterval(range, zero))
        found =false;
    return found;
}

void findZeroRange(const Circuit_F2_1 &ic, int num, Interval range, vector<double> &pos)
{
    double zero;

    if (newtonMid(ic, num, range, zero))
    {
        pos.push_back(zero);
        findZeroRange(ic, num, Interval(range.lower(), zero), pos);
        findZeroRange(ic, num, Interval(zero, range.upper()), pos);
    }
}

void geneChildren(const Circuit_F2_1 &ic, const vector< vector<double> > &pos,
                  size_t num, Circuit_F2_1 child, vector<Circuit_F2_1> &children)
{
    if (num == SIZE_PARM_F2_1)
    {
        children.push_back(child);
        return;
    }

    size_t len = pos[num].size();
    if (len == 0)
    {
        child.set_pi(num, Interval(ic.get_pi(num).lower(), ic.get_pi(num).upper()));
        geneChildren(ic, pos, num+1, child, children);
        return;
    }
    for (size_t i=0; i<len; i++)
    {
        if (i == 0)
        {
            child.set_pi(num, Interval(ic.get_pi(num).lower(), pos[num][i]));
            geneChildren(ic, pos, num+1, child, children);
        } else
        {
            child.set_pi(num, Interval(pos[num][i-1], pos[num][i]));
            geneChildren(ic, pos, num+1, child, children);
        }

        if (i == len-1)
        {
            child.set_pi(num, Interval(pos[num][i], ic.get_pi(num).upper()));
            geneChildren(ic, pos, num+1, child, children);
        }
    }
}

vector<Circuit_F2_1> findZeroF2_1All(const Circuit_F2_1 &ic)
{
    vector< vector<double> > pos(SIZE_PARM_F2_1);
    vector<Circuit_F2_1> children;
    bool found = false;

    for (size_t i=0; i<SIZE_PARM_F2_1; i++)
    {
        pos[i].empty();
        findZeroRange(ic, i, ic.get_pi(i), pos[i]);
        if (pos[i].size() > 0)
            found =true;
        if (pos[i].size() > 1)
            found =true;
    }

    children.empty();
    if (found)
        geneChildren(ic, pos, 0, Circuit_F2_1(), children);

    return children;
}
