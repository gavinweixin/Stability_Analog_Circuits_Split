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

Interval* Interval_7D :: RouthTable(double d) {
    Interval coef[RT_SIZE_F2_1];
    coef[0] = p[1]*p[2]+p[0]*p[2];
    coef[1] = p[5]*p[0]*p[1]*p[2] + p[6]*p[4]*p[1]*p[2] + p[6]*p[0]*p[1]*p[2] - p[6]*p[0]*p[4]*p[3];
    coef[2] = p[6]*p[5]*p[0]*p[4]*p[1]*p[2];

    Interval* RT = new Interval[RT_SIZE_F2_1];
    RT[0] = coef[2]*d*d-coef[1]*d+coef[0];
    RT[1] = coef[1]-coef[2]*2*d;
    RT[2] = coef[2];

    return RT;
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

int Interval_7D::judge()
{
    Interval *RT = RouthTable();
    int positive=0, negtive=0, straddle=0;
    for (int i=0; i<RT_SIZE_F2_1; i++)
    {
        if (RT[i].get_inf()>0) positive++;
        else if (RT[i].get_sup()<0) negtive++;
        else straddle++;
    }
    if (positive==RT_SIZE_F2_1) return 1;
    else if (negtive==RT_SIZE_F2_1) return -1;
    else return 0;
}
