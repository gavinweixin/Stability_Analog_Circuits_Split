//#include "stdafx.h"
#include "Interval.h"
#include <iostream>

Interval :: Interval () : inf(0), sup(0) { }

Interval :: Interval(const double& low, const double& upper) : inf(low), sup(upper) { }

Interval :: Interval (const Interval& orig) : inf(orig.inf), sup(orig.sup) { }

Interval& Interval :: operator = (const Interval& i) {
	inf = i.inf;
	sup = i.sup;
	return *this;
}

Interval& Interval :: operator += (const Interval& i) {
	inf += i.inf;
	sup += i.sup;
	return *this;
}

Interval& Interval :: operator -= (const Interval& i) {
	inf -= i.sup;
	sup -= i.inf;
	return *this;
}

Interval& Interval :: operator *= (const Interval& i) {
	double temp_i(inf);
	double temp_s(sup);
	if ( inf>=0 ) {
		if ( i.inf>=0 ) {
			inf = temp_i*i.inf;
			sup = temp_s*i.sup;
		}
		else if ( i.sup<=0 ) {
			inf = temp_s*i.inf;
			sup = temp_i*i.sup;
		}
		else if ( i.inf<0 && i.sup>=0 ) {
			inf = temp_s*i.inf;
			sup = temp_s*i.sup;
		}
	}
	else if ( sup<=0) {
		if ( i.inf>=0 ) {
			inf = temp_i*i.sup;
			sup = temp_s*i.inf;
		}
		else if ( i.sup<=0 ) {
			inf = temp_s*i.sup;
			sup = temp_i*i.inf;
		}
		else if ( i.inf<0 && i.sup>=0 ) {
			inf = temp_i*i.sup;
			sup = temp_i*i.inf;
		}
	}	
	else if ( inf<0 && sup>=0) {
		if ( i.inf>=0 ) {
			inf = temp_i*i.sup;
			sup = temp_s*i.sup;
		}
		else if ( i.sup<=0 ) {
			inf = temp_s*i.inf;
			sup = temp_i*i.inf;
		}
		else if ( i.inf<0 && i.sup>=0 ) {
			inf = (temp_i*i.sup < temp_s*i.inf) ? temp_i*i.sup : temp_s*i.inf;
			sup = (temp_i*i.inf > temp_s*i.sup) ? temp_i*i.inf : temp_s*i.sup;
		}
	}
	else {
		inf = 0;
		sup = 0;
	}
	return *this;
}


bool Interval :: lt_0 () {
	return (inf>0);
}

bool Interval :: st_0 () {
	return (sup<0);
}

bool Interval :: uncertain() {
	return (inf<=0 && sup>=0);
}

double Interval :: width_cal() {
	return (sup-inf);
}

double Interval :: get_inf() {
	return inf;
}

double Interval :: get_sup() {
	return sup;
}

void Interval :: set_inf(const double& low) {
	this->inf = low;
}

void Interval :: set_sup(const double& upper) {
	this->sup = upper;
}

Interval ::	~Interval() { }

Interval operator + (const Interval& lhs, const Interval& rhs) {
	Interval ret(lhs);
	ret += rhs;
	return ret;
}

Interval operator - (const Interval& lhs, const Interval& rhs) {
	Interval ret(lhs);
	ret -= rhs;
	return ret;
}

Interval operator * (const Interval& lhs, const Interval& rhs) {
	Interval ret(lhs);
	ret *= rhs;
	return ret;
}
