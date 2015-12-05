#include <Circuit_F2_1.h>
#include <cmath>

double newton(Circuit_F2_1 ic, int num, bool isInf, bool &found)
{
    const int MAXITEA = 20;
    vector<Interval> f;
    Interval df;
    double itr;

    ic.set_pi(num, Interval(median(ic.get_pi(num))));
    double eps = median(ic.get_pi(num))/10000;  //set the accuracy
    found = false;
    for (int i=0; i<MAXITEA; i++)
    {
        f = ic.RouthTable();
        df = (ic.Jacobi_cal())[1][num];    //to be modified
        if (isInf)
            itr = median(ic.get_pi(num))-f[1].lower()/df.lower();
        else
            itr = median(ic.get_pi(num))-f[1].upper()/df.upper();
        if (abs(itr-median(ic.get_pi(num)))<eps)
        {
            found = true;
            break;
        }
        ic.set_pi(num, Interval(itr));
    }
    return itr;
}

//bounds of interval not included
bool insideInterval(const Interval &intval, const double &value)
{
    double eps = median(intval)/10000;
//    double eps = width(intval)/10;  //to be modified
    return ((value-intval.lower()>eps) && (intval.upper()-value>eps));
}

vector<double> findZeroF2_1I(const Circuit_F2_1 &ic)
{
    vector<double> pos;
    pos.resize(SIZE_PARM_F2_1);
    bool found;

    for (int i=0; i<SIZE_PARM_F2_1; i++)
    {
        pos[i] = newton(ic, i, true, found);
        if (!found || !insideInterval(ic.get_pi(i), pos[i])/*in(pos[i], ic.get_pi(i))*/)
        {
            pos[i] = newton(ic, i, false, found);
            if (!found || !insideInterval(ic.get_pi(i), pos[i])/*in(pos[i], ic.get_pi(i))*/)
                pos[i] = -1;
        }
    }
    return pos;
}
