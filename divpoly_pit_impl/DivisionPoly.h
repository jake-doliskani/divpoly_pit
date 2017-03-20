#ifndef DIVISIONPOLY_H
#define DIVISIONPOLY_H

#include <NTL/ZZ_pEX.h>
#include <NTL/vector.h>

using namespace NTL;

class DivisionPoly {
public:

    DivisionPoly();
    virtual ~DivisionPoly();

    void compute(Vec<ZZ_pEX> &result,
                 const ZZ_pE &a,
                 const ZZ_pE &b,
                 const ZZ &n,
                 long modulusDegree);

private:

    void computeNextValue(ZZ_pEX &result,
                          const Vec<ZZ_pEX> &auxValuesS,
                          const Vec<ZZ_pEX> &auxValuesT,
                          int currentIndex,
                          int nextIndex);

    void computeFirstValues(Vec<ZZ_pEX> &result,
                            ZZ_pEX &psi2,
                            int msbits);

    void computAuxValues(Vec<ZZ_pEX> &auxValuesS,
                         Vec<ZZ_pEX> &auxValuesT,
                         const Vec<ZZ_pEX> &values);

    long modulusDegree;
    ZZ_pE a;
    ZZ_pE b;
    ZZ_pEX psi2Square;
};

#endif /* DIVISIONPOLY_H */

