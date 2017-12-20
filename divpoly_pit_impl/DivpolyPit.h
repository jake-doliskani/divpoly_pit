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

    bool solvePit(long r);
    void getCycloFactor(ZZ_pEX &factor, long r);

    ZZ_pE a;
    ZZ_pE b;
    
    const double EPSILON = 0.45;
};

#endif /* DIVPOLYPIT_H */

