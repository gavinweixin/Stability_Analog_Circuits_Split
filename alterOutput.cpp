#include <Circuit_F2_1.h>
#include <iostream>
#include <fstream>
#include <iomanip>

bool cubeTest(Circuit cube, size_t pos)
{
    if (pos == cube.SIZE_PARM())
    {
        if (cube.judge())
            return true;
        else
            return false;
    }
    Interval tmp = cube.get_pi(pos);
    bool test;
    cube.set_pi(pos, Interval(tmp.lower()));
    test = cubeTest(cube, pos+1);
    cube.set_pi(pos, Interval(tmp.upper()));
    test = test && cubeTest(cube, pos+1);

    return test;
}

void stableVerify(vector<Circuit> &stable)
{
    bool verified = true;

    for (vector<Circuit>::iterator ivec = stable.begin(); ivec != stable.end(); ++ ivec)
        if (!cubeTest(*ivec, 0))
        {
            verified = false;
            break;
        }
    if (!verified)
        cout << "Fail!" << endl;
    else
        cout << "Pass." <<endl;
}

void normalOut(double totVol, Circuit p)
{
    vector<Circuit> stable;
    vector<Circuit> unstable;
    vector<Circuit> uncertain;

    Judge(totVol,p,stable,unstable,uncertain);
    cout << "#Cube\t" << stable.size() << "\t" << unstable.size() << "\t" << uncertain.size() << endl;

    double vol_stable=0, vol_unstable=0, vol_uncertain=0;
    for (size_t i=0; i<stable.size();    i++) vol_stable    += stable[i].volume_cal();
    for (size_t i=0; i<unstable.size();  i++) vol_unstable  += unstable[i].volume_cal();
    for (size_t i=0; i<uncertain.size(); i++) vol_uncertain += uncertain[i].volume_cal();

    cout << "Vol%\t";
    cout << setprecision(6);
    cout << vol_stable/totVol*100 << "%\t" << vol_unstable/totVol*100 << "%\t"
         << vol_uncertain/totVol*100 << "%" << endl;

    stableVerify(stable);
}

void distribution(double totVol, vector<Circuit> &stable)
{
    vector<size_t> ratio(10, 0);
    double vol = 0.;
    for (vector<Circuit>::iterator ivec = stable.begin(); ivec != stable.end(); ++ ivec)
    {
        ratio[(ivec->volume_cal()/totVol)*10000]++;
        if (ivec->volume_cal()/totVol > 1e-3)
            vol += ivec->volume_cal()/totVol;
    }
    for (size_t i=0; i<10; i++)
        cout << i/100. << "% ~ " << (i+1)/100. << "%:\t" << ratio[i] << endl;
    cout << "sum of vol which is larger than 0.1%: " << vol*100 << "%" << endl;
}

void shiftD(double totVol, Circuit p)
{
    vector<Circuit> stable;
    vector<Circuit> unstable;
    vector<Circuit> uncertain;
    ofstream fout("shift_d");
    double vol_s;
    for (int d=0; d<4000; d+=400)
    {
        stable.clear();
        unstable.clear();
        uncertain.clear();
        Judge(totVol,p,stable,unstable,uncertain,d);
        vol_s = 0;
        for (size_t i=0; i<stable.size(); i++)
            vol_s += stable[i].volume_cal();
        fout << setprecision(6);
        fout << d << "\t" << vol_s/totVol*100 << endl;
    }
    fout.close();
}

void cubePrint(Circuit &p, vector<Circuit> &s, vector<Circuit> &us, vector<Circuit> &uc)
{
    cout << "stable:" << endl;
    for (vector<Circuit>::iterator ivec = s.begin(); ivec != s.end(); ++ ivec)
    {
        for (size_t i = 0; i != p.SIZE_PARM(); ++ i)
            cout << "[" << (*ivec).get_pi(i).lower() << "," << (*ivec).get_pi(i).upper() << "]" << "\t";
        cout << endl;
    }
    cout << "unstable:" << endl;
    for (vector<Circuit>::iterator ivec = us.begin(); ivec != us.end(); ++ ivec)
    {
        for (size_t i = 0; i != p.SIZE_PARM(); ++ i)
           cout << "[" << (*ivec).get_pi(i).lower() << "," << (*ivec).get_pi(i).upper() << "]" << "\t";
        cout << endl;
    }
    cout << "uncertain:" << endl;
    for (vector<Circuit>::iterator ivec = uc.begin(); ivec != uc.end(); ++ ivec)
    {
        for (size_t i = 0; i != p.SIZE_PARM(); ++ i)
            cout << "[" << (*ivec).get_pi(i).lower() << "," << (*ivec).get_pi(i).upper() << "]" << "\t";
        cout << endl;
    }
}
