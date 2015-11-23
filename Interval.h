#ifndef _INTERVAL_H_
#define _INTERVAL_H_

class Interval {
public:
	Interval();
	Interval(const double& low, const double& upper);
	Interval(const Interval& orig);
    Interval& operator = (const Interval& i);
    Interval& operator += (const Interval& i);
    Interval& operator -= (const Interval& i);
    Interval& operator *= (const Interval& i);
	bool lt_0 ();
	bool st_0 ();
	bool uncertain();
	double width_cal();
	double get_inf();
	double get_sup();
	void set_inf(const double& low);
	void set_sup(const double& upper);
	~Interval();
private:
	double inf;
	double sup;
};

Interval operator + (const Interval& lhs, const Interval& rhs);
Interval operator - (const Interval& lhs, const Interval& rhs);
Interval operator * (const Interval& lhs, const Interval& rhs);

#endif
