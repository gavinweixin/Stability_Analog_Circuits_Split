#include "Interval.h"
#include <vector>
using namespace std;

#ifndef _INTERVAL_7D_H_
#define _INTERVAL_7D_H_

class Interval_7D
{
public:
	Interval_7D();
	Interval_7D( const Interval_7D& orig );
	Interval_7D( const vector<Interval>& orig );
	~Interval_7D();
    Interval_7D& operator = (const Interval_7D& i);
	double volume_cal ();
	Interval F();
	vector<Interval> J2_cal();
	Interval& get_pi( int i );
	Interval_7D b_sub_bc ();
private:
	vector<Interval> p;
	double volume;
};

#endif
