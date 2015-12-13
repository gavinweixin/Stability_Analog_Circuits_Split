#include "Circuit_F2_1.h"



Circuit_F2_1 init()
{
    vector<Interval> temp_v;
    Circuit_F2_1 temp;
    temp_v.push_back(Interval(10.00e+3*0.85, 10.00e+3*1.15));
    temp_v.push_back(Interval(10.00e+3*0.85, 10.00e+3*1.15));
    temp_v.push_back(Interval(10.00e+3*0.85, 10.00e+3*1.15));
    temp_v.push_back(Interval(588.1e+3*0.85, 588.1e+3*1.15));
    temp_v.push_back(Interval(245.4   *0.85, 245.4   *1.15));
    temp_v.push_back(Interval(90.87e-9*0.85, 90.87e-9*1.15));
    temp_v.push_back(Interval(90.07e-9*0.85, 90.07e-9*1.15));

    temp = Circuit_F2_1(temp_v);
    return temp;
}

int main()
{
    cout << "Split Method:\t";
    #ifdef SplitMethod_CFBM
        #ifdef findZeroAdded
        cout << "CFBM with findZero" << endl;
        #else
        cout << "CFBM" << endl;
        #endif
    #elif defined SplitMethod_Jacobi
        cout << "Jacobi" << endl;
    #endif

    cout << "Circuit:\t\t";
    #ifdef IC_F2_1
        cout << "F2_1" << endl;
    #endif

    Circuit_F2_1 p = init();
    double totVol = p.volume_cal();
    normalOut(totVol, p);

//    fout.close();
    return 0;
}
