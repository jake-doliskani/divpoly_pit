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
    void cycloFactor(ZZ_pEX &factor, long r);

    ZZ_pE a;
    ZZ_pE b;
    
    const double EPSILON = 0.1;
    const long MIN_THREASHOLD = 17;
};

#endif /* DIVPOLYPIT_H */

