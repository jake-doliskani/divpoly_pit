#include <cstdlib>
#include "DivisionPoly.h"
#include "CycloMod.h"
#include "Util.h"
#include <iostream>
#include <NTL/ZZ_pXFactoring.h>

using namespace std;

void testDivPoly() {
    ZZ_pX f;
    clear(f);
    SetCoeff(f, 2, 1);
    SetCoeff(f, 0, 1);
    ZZ_pE::init(f);
    
    ZZ_pE a, b;
    SetCoeff(a._ZZ_pE__rep, 0, 6);
    SetCoeff(a._ZZ_pE__rep, 1, 6);
    SetCoeff(b._ZZ_pE__rep, 0, 6);
    SetCoeff(b._ZZ_pE__rep, 1, 1);
    cout << f << endl;
    cout << a << endl;
    cout << b << endl;
    cout << "--------------------" << endl;
    Vec<ZZ_pEX> polys;
    polys.SetLength(9);
    ZZ n = to_ZZ(19);
    DivisionPoly divisionPoly;
    divisionPoly.compute(polys, a, b, n, 1000);

    cout << polys[3] << endl;
}

//void testPowerMod() {
//    ZZ_pX f;
//    clear(f);
//    SetCoeff(f, 2, 1);
//    SetCoeff(f, 0, 1);
//    ZZ_pE::init(f);
//    
//    long modulusDegree = 100;
//    
//    ZZ_pEX g, result1, result2;
//    
//    ZZ_pEXModulus F;
//    SetCoeff(g, 0, -1);
//    SetCoeff(g, modulusDegree, 1);
//    build(F, g);
//    
//    random(g, modulusDegree);
//    ZZ e = RandomBits_ZZ(1000);
//    
//    CycloMod cycloMod(modulusDegree);
//    
//    Util util;
//    long start = util.getTimeMillis();
//    cycloMod.powerMod(result1, g, e);
//    cout << util.getTimeMillis() - start << endl;
//    
//    start = util.getTimeMillis();
//    PowerMod(result2, g, e, F);
//    cout << util.getTimeMillis() - start << endl;
//    
//    if (result1 == result2)
//        cout << "OK" << endl;
//    else
//        cout << "Failed" << endl;
//}



int main(int argc, char** argv) {
    ZZ p;
    RandomBits(p, 3);
    NextPrime(p, p);
    ZZ_p::init(p);

    

    return 0;
}

