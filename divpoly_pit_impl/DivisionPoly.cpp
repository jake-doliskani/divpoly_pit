#include "DivisionPoly.h"
#include "CycloMod.h"
#include <iostream>

using namespace std;

DivisionPoly::DivisionPoly() {
}

DivisionPoly::~DivisionPoly() {
    psi2.kill();
}

void DivisionPoly::compute(Vec<ZZ_pEX>& result, 
                           const ZZ_pE& a, 
                           const ZZ_pE& b, 
                           const ZZ& n, 
                           long modulusDegree) {
    
    this->modulusDegree = modulusDegree;
    this->a = a;
    this->b = b;
    
    int msbits = 0;
    long size = NumBits(n);
    msbits |= bit(n, size - 1);
    msbits <<= 1;
    msbits |= bit(n, size - 2);
    msbits <<= 1;
    msbits |= bit(n, size - 3);
    
    if (msbits >= 6)
        msbits = 3;

    computeFirstValues(result, msbits);
    
    Vec<ZZ_pEX> auxValuesS;
    Vec<ZZ_pEX> auxValuesT;

    auxValuesS.SetLength(9);
    auxValuesT.SetLength(9);
    
    int m = 0;
    if (msbits == 3)
        m = size - 3;
    else
        m = size - 4;

    int currentIndex = msbits;
    int nextIndex = 0;
    for (long i = m; i >= 0; i--) {
        currentIndex = currentIndex & 0x3;
        nextIndex = 2 * currentIndex + bit(n, i);
        
        computAuxValues(auxValuesS, auxValuesT, result);
        
        for (int j = -3; j < 6; j++)
            computeNextValue(result[j + 3], auxValuesS, auxValuesT, 
                              currentIndex, nextIndex + j);
        
        currentIndex = nextIndex;
    }
    
    CycloMod cycloMod(modulusDegree);
    
    int odd = bit(n, 0);
    for (int i = 0; i < 9; i++) {
        if ((i % 2 + odd) == 1)
            cycloMod.mulMod(result[i], result[i], psi2);
    }

    auxValuesS.kill();
    auxValuesT.kill();
    psi2SquareEX.kill();
}

void DivisionPoly::evaluate(Vec<ZZ_pE>& result,
                            const ZZ_pE& a,
                            const ZZ_pE& b,
                            const ZZ_pE& x,
                            const ZZ& n) {

    this->a = a;
    this->b = b;
    
    int msbits = 0;
    long size = NumBits(n);
    msbits |= bit(n, size - 1);
    msbits <<= 1;
    msbits |= bit(n, size - 2);
    msbits <<= 1;
    msbits |= bit(n, size - 3);
    
    if (msbits >= 6)
        msbits = 3;

    ZZ_pE psi2;

    computeFirstValues(result, psi2, x, msbits);
    
    Vec<ZZ_pE> auxValuesS;
    Vec<ZZ_pE> auxValuesT;

    auxValuesS.SetLength(9);
    auxValuesT.SetLength(9);
    
    int m = 0;
    if (msbits == 3)
        m = size - 3;
    else
        m = size - 4;

    int currentIndex = msbits;
    int nextIndex = 0;
    for (long i = m; i >= 0; i--) {
        currentIndex = currentIndex & 0x3;
        nextIndex = 2 * currentIndex + bit(n, i);
        
        computAuxValues(auxValuesS, auxValuesT, result);
        
        for (int j = -3; j < 6; j++)
            computeNextValue(result[j + 3], auxValuesS, auxValuesT, 
                              currentIndex, nextIndex + j);
        
        currentIndex = nextIndex;
    }
    
    int odd = bit(n, 0);
    for (int i = 0; i < 9; i++) {
        if ((i % 2 + odd) == 1)
            mul(result[i], result[i], psi2);
    }

    auxValuesS.kill();
    auxValuesT.kill();
}

ZZ_pEX DivisionPoly::getPsi2() const {
    return psi2;
}

void DivisionPoly::computeNextValue(ZZ_pEX &result,
                                    const Vec<ZZ_pEX> &auxValuesS,
                                    const Vec<ZZ_pEX> &auxValuesT,
                                    int currentIndex,
                                    int nextIndex) {
    
    ZZ_pEX temp;
    int j = 0;
    
    CycloMod cycloMod(modulusDegree);
    
    if (nextIndex % 2 == 0) {
        
        j = nextIndex / 2;
        j = j - currentIndex + 3;
        cycloMod.mulMod(temp, auxValuesT[j - 1], auxValuesS[j + 1]);
        cycloMod.mulMod(result, auxValuesT[j + 1], auxValuesS[j - 1]);
        sub(result, temp, result);
        
    } else {
        
        j = (nextIndex - 1) / 2;
        if (j % 2 == 0) {
            
            j = j - currentIndex + 3;
            cycloMod.mulMod(temp, auxValuesT[j], auxValuesS[j + 1]);
            cycloMod.mulMod(result, auxValuesT[j + 1], auxValuesS[j]);
            cycloMod.mulMod(temp, temp, psi2SquareEX);
            sub(result, temp, result);
            
        } else {
            
            j = j - currentIndex + 3;
            cycloMod.mulMod(temp, auxValuesT[j], auxValuesS[j + 1]);
            cycloMod.mulMod(result, auxValuesT[j + 1], auxValuesS[j]);
            cycloMod.mulMod(result, result, psi2SquareEX);
            sub(result, temp, result);
            
        }
    }

    temp.kill();
}

void DivisionPoly::computeNextValue(ZZ_pE &result,
                                    const Vec<ZZ_pE> &auxValuesS,
                                    const Vec<ZZ_pE> &auxValuesT,
                                    int currentIndex,
                                    int nextIndex) {
    
    ZZ_pE temp;
    int j = 0;
    
    if (nextIndex % 2 == 0) {
        
        j = nextIndex / 2;
        j = j - currentIndex + 3;
        mul(temp, auxValuesT[j - 1], auxValuesS[j + 1]);
        mul(result, auxValuesT[j + 1], auxValuesS[j - 1]);
        sub(result, temp, result);
        
    } else {
        
        j = (nextIndex - 1) / 2;
        if (j % 2 == 0) {
            
            j = j - currentIndex + 3;
            mul(temp, auxValuesT[j], auxValuesS[j + 1]);
            mul(result, auxValuesT[j + 1], auxValuesS[j]);
            mul(temp, temp, psi2SquareE);
            sub(result, temp, result);
            
        } else {
            
            j = j - currentIndex + 3;
            mul(temp, auxValuesT[j], auxValuesS[j + 1]);
            mul(result, auxValuesT[j + 1], auxValuesS[j]);
            mul(result, result, psi2SquareE);
            sub(result, temp, result);
            
        }
    }
}


void DivisionPoly::computeFirstValues(Vec<ZZ_pEX>& result, int msbits) {
    
    for (int i = 0; i < 9; i++)
        clear(result[i]);
    
    CycloMod cycloMod(modulusDegree);
    
    // psi2 = 4 * (x^3 + a * x + b)
    SetCoeff(psi2, 3, 1);
    SetCoeff(psi2, 1, a);
    SetCoeff(psi2, 0, b);
    mul(psi2, psi2, 4);
    cycloMod.sqrMod(psi2SquareEX, psi2);
    
    // f1 = 1
    set(result[1]);
    
    // f2 = 1
    set(result[2]);
    
    // 3x^4 + 6Ax^2 + 12Bx − A^2
    SetCoeff(result[3], 4, 3);
    SetCoeff(result[3], 2, 6 * a);
    SetCoeff(result[3], 1, 12 * b);
    SetCoeff(result[3], 0, - a * a);
    cycloMod.reduce(result[3], result[3]);

    // 2*(x^6 + 5Ax^4 + 20Bx^3 − 5A^2x^2 − 4ABx − 8B^2 − A^3)
    SetCoeff(result[4], 6, 1);
    SetCoeff(result[4], 4, 5 * a);
    SetCoeff(result[4], 3, 20 * b);
    SetCoeff(result[4], 2, -5 * a * a);
    SetCoeff(result[4], 1, -4 * a * b);
    SetCoeff(result[4], 0, -8 * b * b - a * a * a);
    mul(result[4], result[4], 2);
    cycloMod.reduce(result[4], result[4]);

    Vec<ZZ_pEX> auxValuesS;
    Vec<ZZ_pEX> auxValuesT;
    auxValuesS.SetLength(9);
    auxValuesT.SetLength(9);
    
    
    computAuxValues(auxValuesS, auxValuesT, result);
    for (int i = 5; i < 9; i++) {
        computeNextValue(result[i], auxValuesS, auxValuesT, 3, i);
        cycloMod.mulMod(auxValuesS[i - 1], result[i - 2], result[i]);
        cycloMod.sqrMod(auxValuesT[i], result[i]);
    }

    if (msbits == 3)
        return;

    if (msbits == 4) {
        
        for (int i = 1; i < 9; i++) {
            result[i - 1] = result[i];
            auxValuesS[i - 1] = auxValuesS[i];
            auxValuesT[i - 1] = auxValuesT[i];
        }
        
        cycloMod.mulMod(auxValuesS[6], result[5], result[7]);
        cycloMod.sqrMod(auxValuesT[7], result[7]);
        computeNextValue(result[8], auxValuesS, auxValuesT, 4, 9);
        return;
    }

    if (msbits == 5) {
        for (int i = 2; i < 9; i++) {
            result[i - 2] = result[i];
            auxValuesS[i - 2] = auxValuesS[i];
            auxValuesT[i - 2] = auxValuesT[i];
        }
        
        cycloMod.mulMod(auxValuesS[5], result[4], result[6]);
        cycloMod.sqrMod(auxValuesT[6], result[6]);
        computeNextValue(result[7], auxValuesS, auxValuesT, 5, 9);
        cycloMod.mulMod(auxValuesS[6], result[5], result[7]);
        cycloMod.sqrMod(auxValuesT[7], result[7]);
        computeNextValue(result[8], auxValuesS, auxValuesT, 5, 10);
        return;
    }
}

void DivisionPoly::computeFirstValues(Vec<ZZ_pE>& result, 
                                      ZZ_pE& psi2,
                                      const ZZ_pE& x,
                                      int msbits) {
    
    for (int i = 0; i < 9; i++)
        clear(result[i]);
    
    // psi2 = 4 * (x^3 + a * x + b)
    psi2 = 4 * (power(x, 3) + a * x + b);
    psi2SquareE = sqr(psi2);
    
    // f1 = 1
    set(result[1]);
    
    // f2 = 1
    set(result[2]);
    
    // 3x^4 + 6ax^2 + 12bx − A^2
    result[3] = 3 * power(x, 4) + 6 * a * sqr(x) + 12 * b * x - sqr(a);

    // 2*(x^6 + 5ax^4 + 20bx^3 − 5a^2x^2 − 4abx − 8b^2 − a^3)
    result[4] = power(x, 6) + 5 * a * power(x, 4) + 20 * b * power(x, 3) -
            5 * sqr(a * x) - 4 * a * b * x - 8 * sqr(b) - power(a, 3);
    result[4] *= 2;

    Vec<ZZ_pE> auxValuesS;
    Vec<ZZ_pE> auxValuesT;
    auxValuesS.SetLength(9);
    auxValuesT.SetLength(9);
    
    
    computAuxValues(auxValuesS, auxValuesT, result);
    for (int i = 5; i < 9; i++) {
        computeNextValue(result[i], auxValuesS, auxValuesT, 3, i);
        mul(auxValuesS[i - 1], result[i - 2], result[i]);
        sqr(auxValuesT[i], result[i]);
    }

    if (msbits == 3)
        return;

    if (msbits == 4) {
        
        for (int i = 1; i < 9; i++) {
            result[i - 1] = result[i];
            auxValuesS[i - 1] = auxValuesS[i];
            auxValuesT[i - 1] = auxValuesT[i];
        }
        
        mul(auxValuesS[6], result[5], result[7]);
        sqr(auxValuesT[7], result[7]);
        computeNextValue(result[8], auxValuesS, auxValuesT, 4, 9);
        return;
    }

    if (msbits == 5) {
        for (int i = 2; i < 9; i++) {
            result[i - 2] = result[i];
            auxValuesS[i - 2] = auxValuesS[i];
            auxValuesT[i - 2] = auxValuesT[i];
        }
        
        mul(auxValuesS[5], result[4], result[6]);
        sqr(auxValuesT[6], result[6]);
        computeNextValue(result[7], auxValuesS, auxValuesT, 5, 9);
        mul(auxValuesS[6], result[5], result[7]);
        sqr(auxValuesT[7], result[7]);
        computeNextValue(result[8], auxValuesS, auxValuesT, 5, 10);
        return;
    }
}

void DivisionPoly::computAuxValues(Vec<ZZ_pEX>& auxValuesS, 
                                   Vec<ZZ_pEX>& auxValuesT, 
                                   const Vec<ZZ_pEX>& values) {
    CycloMod cycloMod(modulusDegree);
    
    for (int i = 1; i < values.length() - 1; i++) {
        cycloMod.mulMod(auxValuesS[i], values[i - 1], values[i + 1]);
        cycloMod.sqrMod(auxValuesT[i], values[i]);
    }
}

void DivisionPoly::computAuxValues(Vec<ZZ_pE>& auxValuesS, 
                                   Vec<ZZ_pE>& auxValuesT, 
                                   const Vec<ZZ_pE>& values) {
    for (int i = 1; i < values.length() - 1; i++) {
        mul(auxValuesS[i], values[i - 1], values[i + 1]);
        sqr(auxValuesT[i], values[i]);
    }
}