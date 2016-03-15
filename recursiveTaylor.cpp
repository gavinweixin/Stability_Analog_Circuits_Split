#include "Circuit_F2_1.h"

Interval bound(ex fun, const vector<symbol> &sym, const vector<Interval> &param)
{
    size_t n = sym.size();

    bool constant = true;
    for (size_t i=0; i<n; i++)
        if (fun.has(sym[i]))
        {
            constant = false;
            break;
        }
    if (constant)
        return Interval(to_double(ex_to<numeric>(fun)));

    exmap p0;
    vector<double> p1;
    for (size_t i=0; i<n; i++)
    {
        p0[sym[i]] = (param[i].lower()+param[i].upper())/2;
        p1.push_back((param[i].upper()-param[i].lower())/2);
    }

    Interval b(to_double(ex_to<numeric>(subs(fun, p0))));
    for (size_t i=0; i<n; i++)
    {
        b = b + p1[i]*to_double(ex_to<numeric>(subs(fun.diff(sym[i]), p0)))*Interval(-1,1);
        b = b + 0.5*p1[i]*p1[i]*Interval(0,1)*bound((fun.diff(sym[i])).diff(sym[i]), sym, param);
    }
    for (size_t i=0; i<n-1; i++)
        for (size_t j=i+1; j<n; j++)
            b = b + p1[i]*p1[j]*Interval(-1,1)*bound((fun.diff(sym[i])).diff(sym[j]), sym, param);

    return b;
}
