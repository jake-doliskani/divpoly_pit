#include "CycloMod.h"

CycloMod::CycloMod(long modulusDegree) {
    this->modulusDegree = modulusDegree;
}

CycloMod::~CycloMod() {
}

void CycloMod::mulMod(ZZ_pEX& result, 
                          const ZZ_pEX& f, 
                          const ZZ_pEX& g) {
    
    ZZ_pE temp1;
    ZZ_pE temp2;
    ZZ_pEX tempResult;
    clear(tempResult);
    
    for (long k = 0; k < modulusDegree; k++) {
        temp2 = 0;
        for (long i = 0; i <= k; i++) {
            mul(temp1, coeff(f, i), coeff(g, k - i));
            add(temp2, temp2, temp1);
        }

        for (long i = k + 1; i < modulusDegree; i++) {
            mul(temp1, coeff(f, i), coeff(g, modulusDegree + k - i));
            add(temp2, temp2, temp1);
        }
        
        SetCoeff(tempResult, k, temp2);
    }
    
    result = tempResult;
    tempResult.kill();

//    mul(result, f, g);
//    reduce(result, result);
}


void CycloMod::invMod(ZZ_pEX& result, const ZZ_pEX& f) {
    ZZ_pEX modulus;
    clear(modulus);
    SetCoeff(modulus, 0, -1);
    SetCoeff(modulus, modulusDegree, 1);
    
    InvMod(result, f, modulus);
    
    modulus.kill();
}


void CycloMod::sqrMod(ZZ_pEX& result, 
                          const ZZ_pEX& f) {
    sqr(result, f);
    reduce(result, result);
}

void CycloMod::reduce(ZZ_pEX& result, const ZZ_pEX& a) {
    if (deg(a) < modulusDegree) {
        result = a;
        return;
    }
    
    long degree = deg(a);
    long j = 0;
    ZZ_pEX temp;
    clear(temp);

    for (long i = 0; i <= degree; i++) {
        SetCoeff(temp, j, coeff(a, i) + coeff(temp, j));
        j++;
        
        if (j == modulusDegree)
            j = 0;
    }
    
    result = temp;
    temp.kill();
}

