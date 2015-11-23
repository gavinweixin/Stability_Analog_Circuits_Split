//#include "stdafx.h"
#include "Interval_7D.h"
#include <iostream>

Interval_7D :: Interval_7D() : volume(0) { }

Interval_7D :: Interval_7D( const Interval_7D &orig) : p(orig.p), volume(orig.volume) { }

Interval_7D :: Interval_7D( const vector<Interval>& orig ) : p(orig), volume(0) { }

Interval_7D :: ~Interval_7D() { }

Interval_7D& Interval_7D :: operator = (const Interval_7D& i) {
	p = i.p;
	volume = i.volume;

	return *this;
}

double Interval_7D :: volume_cal() {
	volume = p[0].width_cal() * p[1].width_cal() * p[2].width_cal() * p[3].width_cal() * p[4].width_cal() * p[5].width_cal() * p[6].width_cal();
	return volume;
}

Interval Interval_7D :: F() {
//	cout << "[" << (p[1]*p[2]).get_inf() << "," << (p[1]*p[2]).get_sup() << "]" << ",";
//	cout << "[" << (p[0]*p[5]).get_inf() << "," << (p[0]*p[5]).get_sup() << "]" << ",";
//	cout << "[" << (p[4]*p[6]).get_inf() << "," << (p[4]*p[6]).get_sup() << "]" << ",";
//	cout << "[" << (p[0]*p[6]).get_inf() << "," << (p[0]*p[6]).get_sup() << "]" << ",";
//	cout << "[" << (p[3]*p[4]).get_inf() << "," << (p[3]*p[4]).get_sup() << "]" << ",";
//	cout << endl;
//	cout << "[" << (p[0]*p[1]*p[2]*p[5]).get_inf() << "," << (p[0]*p[1]*p[2]*p[5]).get_sup() << "]" << ",";
//	cout << "[" << (p[1]*p[2]*p[4]*p[6]).get_inf() << "," << (p[1]*p[2]*p[4]*p[6]).get_sup() << "]" << ",";
//	cout << "[" << (p[0]*p[1]*p[2]*p[6]).get_inf() << "," << (p[0]*p[1]*p[2]*p[6]).get_sup() << "]" << ",";
//	cout << "[" << (p[0]*p[3]*p[4]*p[6]).get_inf() << "," << (p[0]*p[3]*p[4]*p[6]).get_sup() << "]" << ",";
//	cout << endl;
//	cout << "[" << (p[0]*p[1]*p[2]*p[5] + p[1]*p[2]*p[4]*p[6]).get_inf() << "," << (p[0]*p[1]*p[2]*p[5] + p[1]*p[2]*p[4]*p[6]).get_sup() << "]" << ",";
//	cout << "[" << (p[0]*p[1]*p[2]*p[6] - p[0]*p[3]*p[4]*p[6]).get_inf() << "," << (p[0]*p[1]*p[2]*p[6] - p[0]*p[3]*p[4]*p[6]).get_sup() << "]" << ",";
//	cout << endl;
	return (p[0]*p[1]*p[2]*p[5] + p[1]*p[2]*p[4]*p[6] + p[0]*p[1]*p[2]*p[6] - p[0]*p[3]*p[4]*p[6]);
}

vector<Interval> Interval_7D :: J2_cal() {
	vector<Interval> temp;
	Interval temp_i;
	temp_i = p[5]*p[1]*p[2] + p[6]*p[1]*p[2] - p[6]*p[4]*p[3];
	temp.push_back(temp_i);
    temp_i = p[5]*p[0]*p[2] + p[6]*p[0]*p[2] + p[6]*p[4]*p[2];
	temp.push_back(temp_i);
	temp_i = p[5]*p[0]*p[1] + p[6]*p[0]*p[1] + p[6]*p[4]*p[1];
	temp.push_back(temp_i);
    temp_i = Interval(0,0) - p[6]*p[0]*p[4];
//	temp_i.set_inf(-temp_i.get_inf());
//	temp_i.set_sup(-temp_i.get_sup());
	temp.push_back(temp_i);
	temp_i = p[6]*p[1]*p[2] - p[6]*p[0]*p[3];
	temp.push_back(temp_i);
	temp_i = p[0]*p[1]*p[2];
	temp.push_back(temp_i);
	temp_i = p[0]*p[1]*p[2] - p[0]*p[4]*p[3] + p[4]*p[1]*p[2];
	temp.push_back(temp_i);

	return temp;
}


Interval& Interval_7D :: get_pi( int i ) {
	return p[i];
}

Interval_7D Interval_7D :: b_sub_bc() {
	vector<double> inf(7);
	vector<double> sup(7);
	vector<Interval> temp_v;
	Interval_7D temp;

	for ( int i = 0; i != 7; ++ i ) {
        inf[i] = p[i].get_inf() - ( p[i].get_sup()+p[i].get_inf() )/2;
        sup[i] = p[i].get_sup() - ( p[i].get_sup()+p[i].get_inf() )/2;
		temp_v.push_back(Interval(inf[i],sup[i]));
	}
	
	temp = Interval_7D(temp_v);

	return temp;
}
