#include <cstdlib>
#include "DivisionPoly.h"
#include "CycloMod.h"
#include "Util.h"
#include "DivpolyPit.h"
#include <iostream>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZX.h>

using namespace std;

void testPIT() {
    ZZ_pX f;
    clear(f);
    SetCoeff(f, 2, 1);
    SetCoeff(f, 0, 1);
    ZZ_pE::init(f);
    
    ZZ_pE a, b;
    SetCoeff(a._ZZ_pE__rep, 0, 1);
    SetCoeff(a._ZZ_pE__rep, 1, 0);
    SetCoeff(b._ZZ_pE__rep, 0, 0);
    SetCoeff(b._ZZ_pE__rep, 1, 0);
    
    cout << f << endl;
    cout << a << endl;
    cout << b << endl;
    cout << "--------------------" << endl;
    
    DivpolyPit divpolyPit;
    Util util;
    long start = util.getTimeMillis();
    bool ss = divpolyPit.isSuperSingular(a, b);
    cout << util.getTimeMillis() - start << endl;
    
    if (ss)
        cout << "supersingular" << endl;
    else
        cout << "ordinary" << endl;
}

int main(int argc, char** argv) {
    ZZ p = to_ZZ("2754956713021720098890780431368421677372787959271790387477"
            "317890624118510598679299103597525957767387725819347045113854046"
            "173588366530665763998657581236594775145520555464958540941548783"
            "443208819825034730307049044197015309642918491073582811802566074"
            "653557629404118309436594332574045232220927558523596423659928957"
            "87403940341806275640259285102630974616178256580182367");
    
    ZZ_p::init(p);

    testPIT();

    return 0;
}

