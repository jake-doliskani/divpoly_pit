#include "DivpolyPit.h"
#include "Util.h"
#include "DivisionPoly.h"
#include "CycloMod.h"
#include <iostream>
#include <NTL/ZZ_pXFactoring.h>


using namespace std;


DivpolyPit::DivpolyPit() {
}

DivpolyPit::~DivpolyPit() {
}

bool DivpolyPit::isSupersingular(const ZZ_pE& a, const ZZ_pE& b) {
    this->a = a;
    this->b = b;
    
    double epsilon = 0.1;
    
    long d = 0;
    long r = 3;
    ZZ p2 = sqr(ZZ_p::modulus());
    
    Util util;
    while (true) {
        long a = rem(p2, r);
        d = util.findOrder(a, r);
        
        if (pow(d, 2 - epsilon) > r)
            break;
        
        r = NextPrime(r + 1);
    }
    
    return solvePit(3, d);
    
}

bool DivpolyPit::solvePit(long r, long d) {
    Vec<ZZ_pE> values;
    values.SetLength(9);
    
    DivisionPoly divisionPoly;
    ZZ_pE constantTerm;
    constantTerm = 0;
    divisionPoly.evaluate(values, a, b, constantTerm, ZZ_p::modulus());
    constantTerm = values[3];
    
    if (constantTerm != 1 && constantTerm != -1)
        return false;
    
    divisionPoly.evaluate(values, a, b, random_ZZ_pE(), ZZ_p::modulus());
    
    if (values[3] != constantTerm)
        return false;
    
    Vec<ZZ_pEX> polys;
    polys.SetLength(9);
    divisionPoly.compute(polys, a, b, ZZ_p::modulus(), r);
    
    if (deg(polys[3]) > 0)
        return false;
    
    if (coeff(polys[3], 0) != constantTerm)
        return false;
    
    ZZ_pEX xp2;
    ZZ_pEX temp;
    ZZ_pEX psi2;
    
    SetX(temp);
    clear(xp2);
    SetCoeff(xp2, rem(sqr(ZZ_p::modulus()), r), 1);
    psi2 = divisionPoly.getPsi2();

    // test \pis_{p - 1}\psi_{p + 1} = x - x^{p^2}
    sub(xp2, temp, xp2);
    CycloMod cycloMod(r);
    cycloMod.mulMod(temp, polys[2], polys[4]);
    cycloMod.invMod(psi2, psi2);
    cycloMod.mulMod(temp, temp, psi2);
    if (temp != xp2)
        return false;
    
    return true;
}

void DivpolyPit::cycloFactor(ZZ_pEX& factor, long r, long d) {
    ZZ_pX cycloPoly;
    for (long i = 0; i < r; i++)
        SetCoeff(cycloPoly, i, 1);
    
    // compute x^p mod cycloPoly
    long t = rem(ZZ_p::modulus(), r);
    ZZ_pX temp;
    temp = 0;
    SetCoeff(temp, t, 1);
    rem(temp, temp, cycloPoly);
    
    Vec<ZZ_pX> factors;
    EDF(factors, cycloPoly, temp, d);
    
    conv(factor, factors[0]);
    
    cycloPoly.kill();
    temp.kill();
    factors.kill();
}