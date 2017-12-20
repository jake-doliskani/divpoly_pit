#include "DivpolyPit.h"
#include "Util.h"
#include "DivisionPoly.h"
#include "CycloMod.h"
#include <iostream>
#include <NTL/ZZ_pEXFactoring.h>


using namespace std;

DivpolyPit::DivpolyPit() {
}

DivpolyPit::~DivpolyPit() {
}

bool DivpolyPit::isSupersingular(const ZZ_pE& a, const ZZ_pE& b) {
    this->a = a;
    this->b = b;

    long r = NumBits(ZZ_p::modulus());
    r = r / 3;
    r = NextPrime(r);

    Util util;
    while (true) {
        long d = rem(ZZ_p::modulus(), r);
        d = util.findOrder(d, r);

        if (d == r - 1 && pow(d, 2 - EPSILON) > r) {
            break;
        }

        r = NextPrime(r + 1);
    }

    return solvePit(r);
}

bool DivpolyPit::solvePit(long r) {
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
    
    ZZ_pEX cycloFactor;
    getCycloFactor(cycloFactor, r);
    rem(temp, temp, cycloFactor);
    rem(xp2, xp2, cycloFactor);
    
    if (temp != xp2)
        return false;

    return true;
}

void DivpolyPit::getCycloFactor(ZZ_pEX& factor, long r) {
    ZZ_pEX cycloPoly;
    for (long i = 0; i < r; i++)
        SetCoeff(cycloPoly, i, 1);

    // compute x^p mod cycloPoly
    long t = rem(ZZ_p::modulus(), r);
    ZZ_pEX temp;
    temp = 0;
    SetCoeff(temp, t, 1);
    rem(temp, temp, cycloPoly);

    Vec<ZZ_pEX> factors;
    EDF(factors, cycloPoly, temp, (r - 1) / 2);

    factor = factors[0];

    cycloPoly.kill();
    temp.kill();
    factors.kill();
}