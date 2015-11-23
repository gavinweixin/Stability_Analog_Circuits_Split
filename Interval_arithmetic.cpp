// Interval arithmetic.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <vector>
//#include <cmath>
#include "Interval_7D.h"
using namespace std;

const double epsilon = 0.0001;
double V;
int counter = 0;

inline double abs(double i) {
	return i > 0 ? i : -i;
}

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

pair<Interval_7D, Interval_7D> CFBM ( Interval_7D& orig ) {
    //RouthTable here has only one entry
    double tmp, min_v=2*1;
    int min_i=0;
    double wid_RouthT = orig.F().width_cal();
    pair<Interval_7D, Interval_7D> children;

    for ( int i = 0; i != 7; ++ i )
    {
        children = Bisect_j(orig,i);
        Interval RT = children.first.F();
        tmp = RT.width_cal()/wid_RouthT;
        RT = children.second.F();
        tmp += RT.width_cal()/wid_RouthT;
        if (tmp <= min_v)
        {
            min_v = tmp;
            min_i = i;
        }
    }
    counter++;
    return Bisect_j(orig,min_i);
}

void Judge ( Interval_7D& parent, vector<Interval_7D>& s, vector<Interval_7D>& us, vector<Interval_7D>& uc ) {
    pair<Interval_7D, Interval_7D> children;
	if ( parent.F().lt_0() ) {
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
        children = CFBM( parent );
		Judge(children.first,s,us,uc);
		Judge(children.second,s,us,uc);
	}
}

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
    V = p.volume_cal();
    cout << V << endl;
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
    cout << vol_stable/V << "\t" << vol_unstable/V << "\t" << vol_uncertain/V << endl;
	return 0;
}

