// Interval arithmetic.cpp : Defines the entry point for the console application.
//

//#include "stdafx.h"
#include <iostream>
#include <vector>
#include <cmath>
#include "Interval_7D.h"
using namespace std;

const double epsilon = 0.0001;
const double V = 274.9464231;
int counter = 0;

/*inline double abs(double i) {
	return i > 0 ? i : -i;
}*/

inline double max(double m1, double m2) {
	return m1 > m2 ? m1 : m2;
}

pair<Interval_7D, Interval_7D> Bisect_j ( Interval_7D& orig, int j ) {
	Interval_7D temp1(orig), temp2(orig);
	double mid_j;
	mid_j = ( orig.get_pi(j).get_inf() + orig.get_pi(j).get_sup() ) / 2;
	temp1.get_pi(j).set_sup(mid_j);
	temp2.get_pi(j).set_inf(mid_j);
	pair<Interval_7D, Interval_7D> temp;
	temp = make_pair(temp1, temp2);
	return temp;
}

pair<Interval_7D, Interval_7D> Jacobi ( Interval_7D& orig ) {
	double max_df = 0;
    int max_j = 0;
	double d;
	vector<Interval> df2;
	Interval_7D b_subt_bc = orig.b_sub_bc();
	df2 = orig.J2_cal();
	for ( int i = 0; i != 7; ++ i ) {
		df2[i] = df2[i] * b_subt_bc.get_pi(i);
		d = max( abs(df2[i].get_inf()), abs(df2[i].get_sup()) );
		if ( d > max_df ) {
			max_df = d;
			max_j = i;
		}
	}
	counter++;
	return Bisect_j(orig,max_j);
}

void Judge ( Interval_7D& parent, vector<Interval_7D>& s, vector<Interval_7D>& us, vector<Interval_7D>& uc ) {
    pair<Interval_7D, Interval_7D> children;
	if ( parent.F().lt_0() ) {
//		cout << parent.F().get_inf() << endl;
		s.push_back( parent );
		return;
	}
	else if ( parent.F().st_0() ) {
		us.push_back( parent );
		return;
	}
	else if ( parent.volume_cal() / V < epsilon ) {
		uc.push_back( parent );
		return;
	}
	else {
        children = Jacobi( parent );
		Judge(children.first,s,us,uc);
		Judge(children.second,s,us,uc);
	}
}

/*
the Jacobi Matrix:
J = 
[                              R3,                              R3,                         R1 + R2,          0,                   0,               0,                                0 ]
[ C1*R2*R3 + C2*R2*R3 - C2*R13*R4, C1*R1*R3 + C2*R1*R3 + C2*R13*R3, C1*R1*R2 + C2*R1*R2 + C2*R13*R2, -C2*R1*R13, C2*R2*R3 - C2*R1*R4,        R1*R2*R3, R1*R2*R3 - R1*R13*R4 + R13*R2*R3 ]
[                 C1*C2*R13*R2*R3,                 C1*C2*R1*R13*R3,                 C1*C2*R1*R13*R2,          0,      C1*C2*R1*R2*R3, C2*R1*R13*R2*R3,                  C1*R1*R13*R2*R3 ]
in the term of pi
J =
[                                             p[2],                                             p[2],                                      p[0] + p[1],               0,                               0,              0,                                                0 ]
[ p[5]*p[1]*p[2] + p[6]*p[1]*p[2] - p[6]*p[4]*p[3], p[5]*p[0]*p[2] + p[6]*p[0]*p[2] + p[5]*p[4]*p[2], p[5]*p[0]*p[1] + p[6]*p[0]*p[1] + p[6]*p[4]*p[1], -p[6]*p[0]*p[4], p[6]*p[1]*p[2] - p[6]*p[0]*p[3], p[0]*p[1]*p[2], p[0]*p[1]*p[2] - p[0]*p[4]*p[3] + p[4]*p[1]*p[2] ] 
[                         p[5]*p[6]*p[4]*p[1]*p[2],                         p[5]*p[6]*p[0]*p[4]*p[2],                         p[5]*p[6]*p[0]*p[4]*p[1],               0,        p[5]*p[6]*p[0]*p[1]*p[2], p[6]*p[0]*p[4]*p[1]*p2],                p[5]*p[0]*p[4]*p[1]*p[2] ]
*/

Interval_7D init() {
	vector<Interval> temp_v;
	Interval_7D temp;
    temp_v.push_back(Interval(10.00e+3*0.85, 10.00e+3*1.15));
    temp_v.push_back(Interval(10.00e+3*0.85, 10.00e+3*1.15));
    temp_v.push_back(Interval(10.00e+3*0.85, 10.00e+3*1.15));
    temp_v.push_back(Interval(588.1e+3*0.85, 588.1e+3*1.15));
    temp_v.push_back(Interval(245.4   *0.85, 245.4   *1.15));
    temp_v.push_back(Interval(90.87e-9*0.85, 90.87e-9*1.15));
    temp_v.push_back(Interval(90.07e-9*0.85, 90.07e-9*1.15));
//	temp_v.push_back(Interval(2.498,3.765));
//	temp_v.push_back(Interval(-1.874,6.984));
//	temp_v.push_back(Interval(-6.831,-2.349));
//	temp_v.push_back(Interval(5.874,6.294));
//	temp_v.push_back(Interval(1.932,3.132));
//	temp_v.push_back(Interval(-3.132,-1.132));
//	temp_v.push_back(Interval(7.018,7.918));
	/*temp_v.push_back(Interval(5.018,6.932));
	temp_v.push_back(Interval(1.918,7.232));
	temp_v.push_back(Interval(2.438,7.912));
	temp_v.push_back(Interval(9.132,9.794));
	temp_v.push_back(Interval(2.782,5.134));
	temp_v.push_back(Interval(4.362,6.104));
	temp_v.push_back(Interval(3.302,8.003));*/
	/*temp_v.push_back(Interval(0.756,1.064));
	temp_v.push_back(Interval(6.026,8.964));
	temp_v.push_back(Interval(7.936,8.014));
	temp_v.push_back(Interval(2.036,4.384));
	temp_v.push_back(Interval(0.746,3.913));
	temp_v.push_back(Interval(3.906,7.025));
	temp_v.push_back(Interval(7.606,7.918));*/
	temp = Interval_7D(temp_v);
	return temp;
}

int main()
{
	vector<Interval_7D> stable;
	vector<Interval_7D> unstable;
	vector<Interval_7D> uncertain;
	Interval_7D p;
	p = init();
	double temp;
	temp = p.volume_cal();
	cout << temp << endl;
	Judge(p,stable,unstable,uncertain);
//	cout << "stable:" << endl;
//	for ( vector<Interval_7D>::iterator ivec = stable.begin(); ivec != stable.end(); ++ ivec ) {
//		for ( int i = 0; i != 7; ++ i )
//			cout << "[" << (*ivec).get_pi(i).get_inf() << "," << (*ivec).get_pi(i).get_sup() << "]" << "\t";
//		cout << endl;
//	}
//	cout << "unstable:" << endl;
//	for ( vector<Interval_7D>::iterator ivec = unstable.begin(); ivec != unstable.end(); ++ ivec ) {
//		for ( int i = 0; i != 7; ++ i )
//			cout << "[" << (*ivec).get_pi(i).get_inf() << "," << (*ivec).get_pi(i).get_sup() << "]" << "\t";
//		cout << endl;
//	}
//	cout << "uncertain:" << endl;
// 	for ( vector<Interval_7D>::iterator ivec = uncertain.begin(); ivec != uncertain.end(); ++ ivec ) {
//		for ( int i = 0; i != 7; ++ i )
//			cout << "[" << (*ivec).get_pi(i).get_inf() << "," << (*ivec).get_pi(i).get_sup() << "]" << "\t";
//		cout << endl;
//	}
	cout << counter << endl;
	cout << stable.size() << "\t" << unstable.size() << "\t" << uncertain.size() << endl;
    double vol_stable=0, vol_unstable=0, vol_uncertain=0;
    for (size_t i=0; i<stable.size(); i++) vol_stable += stable[i].volume_cal();
    for (size_t i=0; i<unstable.size(); i++) vol_unstable += unstable[i].volume_cal();
    for (size_t i=0; i<uncertain.size(); i++) vol_uncertain += uncertain[i].volume_cal();
    cout << vol_stable << "\t" << vol_unstable << "\t" << vol_uncertain << endl;
    cout << vol_stable/temp << "\t" << vol_unstable/temp << "\t" << vol_uncertain/temp << endl;
	return 0;
}

