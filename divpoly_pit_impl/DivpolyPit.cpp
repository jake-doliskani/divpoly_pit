#include "DivpolyPit.h"

DivpolyPit::DivpolyPit() {
}

DivpolyPit::~DivpolyPit() {
}

bool DivpolyPit::isSuperSingular(const ZZ_pE& a, const ZZ_pE& b) {
    this->a = a;
    this->b = b;
    
    double eplsilon = 0.1;
    
    long d = 0;
    long r = 5;
    ZZ p = ZZ_p::modulus();
    
    
    
    
}

void DivpolyPit::solvePit() {
    
}