#ifndef _INTERVAL_H_
#define _INTERVAL_H_

class Interval
{
public:
    Interval();
    Interval(const double& value);
    Interval(const double& low, const double& upper);
    Interval(const Interval& orig);
    Interval& operator = (const Interval& i);
    Interval& operator += (const Interval& i);
    Interval& operator -= (const Interval& i);
    Interval& operator *= (const Interval& i);
    bool lt_0 () const;
    bool st_0 () const;
//    bool uncertain() const;
    double width_cal() const;
    double get_inf() const;
    double get_sup() const;
    double get_mid() const;
    void set_inf(const double& low);
    void set_sup(const double& upper);
    ~Interval();
private:
    double inf;
    double sup;
};

Interval operator + (const Interval& lhs, const Interval& rhs);
Interval operator + (const Interval& lhs, const double& rhs);
Interval operator - (const Interval& lhs, const Interval& rhs);
Interval operator - (const Interval& lhs, const double& rhs);
Interval operator * (const Interval& lhs, const Interval& rhs);
Interval operator * (const Interval& lhs, const double& rhs);

#endif
