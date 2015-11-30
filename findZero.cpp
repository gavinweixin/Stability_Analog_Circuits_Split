#include <Interval_7D.h>
#include <cmath>

double newton(Interval_7D ic, int num, bool isInf, bool &found)
{
    const int MAXITEA = 20;
    vector<Interval> f;
    Interval df;
    double itr;

    ic.set_pi(num, Interval(ic.get_pi(num).get_mid()));
    double eps = ic.get_pi(num).get_mid()/10000;  //set the accuracy
    found = false;
    for (int i=0; i<MAXITEA; i++)
    {
        f = ic.RouthTable();
        df = (ic.Jacobi_cal())[1][num];    //to be modified
        if (isInf)
            itr = ic.get_pi(num).get_mid()-f[1].get_inf()/df.get_inf();
        else
            itr = ic.get_pi(num).get_mid()-f[1].get_sup()/df.get_sup();
        if (abs(itr-ic.get_pi(num).get_mid())<eps)
        {
            found = true;
            break;
        }
        ic.set_pi(num, Interval(itr));
    }
    return itr;
}

//terminals of interval not included
bool insideInterval(const Interval &intval, const double &value)
{
    double eps = intval.get_mid()/10000;
    return ((value-intval.get_inf()>eps) && (intval.get_sup()-value>eps));
}

vector<double> findZeroF2_1I(const Interval_7D &ic)
{
    vector<double> pos;
    pos.resize(SIZE_PARM_F2_1);
    bool found;

    for (int i=0; i<SIZE_PARM_F2_1; i++)
    {
        pos[i] = newton(ic, i, true, found);
        if (!found || !insideInterval(ic.get_pi(i), pos[i]))
        {
            pos[i] = newton(ic, i, false, found);
            if (!found || !insideInterval(ic.get_pi(i), pos[i])) pos[i] = -1;
        }
    }
    return pos;
}

