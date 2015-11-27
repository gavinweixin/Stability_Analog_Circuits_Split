#include <Interval_7D.h>
#include <cmath>

const double eps = 1e-12;

double newton(Interval_7D ic, int num, bool isInf, bool &found)
{
    const int MAXITEA = 20;
    Interval *f, df;
    double itr;

    ic.set_pi(num, Interval(ic.get_pi(num).get_mid()));
//    eps = ic.get_pi(num).get_mid()/100;
    found = false;
    for (int i=0; i<MAXITEA; i++)
    {
        f = ic.RouthTable();
        df = (ic.J2_cal())[num];
        if (isInf)
            itr = ic.get_pi(num).get_mid()-f[1].get_inf()/df.get_inf();
        else
            itr = ic.get_pi(num).get_mid()-f[1].get_sup()/df.get_sup();
        delete[] f;
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
bool insideInterval(const Interval &intval, const int &value)
{
    return ((value-intval.get_inf()>eps) && (intval.get_sup()-value>eps));
}

double* findZeroF2_1I(const Interval_7D &ic)
{
    double *pos = new double[SIZE_PARM_F2_1];
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

