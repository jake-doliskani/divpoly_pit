#ifndef DIVPOLYPIT_H
#define DIVPOLYPIT_H

#include <NTL/ZZ_pEX.h>


using namespace NTL;

class DivpolyPit {
    
public:
    
    DivpolyPit();
    virtual ~DivpolyPit();
    
    bool isSuperSingular(const ZZ_pE &a, const ZZ_pE &b);
    void solvePit();
    
private:

    ZZ_pE a;
    ZZ_pE b;
};

#endif /* DIVPOLYPIT_H */

