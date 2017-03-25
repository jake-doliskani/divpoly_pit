#ifndef DIVPOLYPIT_H
#define DIVPOLYPIT_H

#include <NTL/ZZ_pEX.h>


using namespace NTL;

class DivpolyPit {
    
public:
    
    DivpolyPit();
    virtual ~DivpolyPit();
    
    bool isSupersingular(const ZZ_pE &a, const ZZ_pE &b);
    
private:

    bool solvePit(long r, long d);
    void cycloFactor(ZZ_pEX &factor, long r, long d);

    ZZ_pE a;
    ZZ_pE b;
};

#endif /* DIVPOLYPIT_H */

