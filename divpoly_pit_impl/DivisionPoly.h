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

    void evaluate(Vec<ZZ_pE> &result,
                  const ZZ_pE &a,
                  const ZZ_pE &b,
                  const ZZ_pE &x,
                  const ZZ &n);

    ZZ_pEX getPsi2() const;

private:

    void computeNextValue(ZZ_pEX &result,
                          const Vec<ZZ_pEX> &auxValuesS,
                          const Vec<ZZ_pEX> &auxValuesT,
                          int currentIndex,
                          int nextIndex);

    void computeNextValue(ZZ_pE &result,
                          const Vec<ZZ_pE> &auxValuesS,
                          const Vec<ZZ_pE> &auxValuesT,
                          int currentIndex,
                          int nextIndex);

    void computeFirstValues(Vec<ZZ_pEX> &result, int msbits);

    void computeFirstValues(Vec<ZZ_pE> &result,
                            ZZ_pE &psi2,
                            const ZZ_pE& x,
                            int msbits);

    void computAuxValues(Vec<ZZ_pEX> &auxValuesS,
                         Vec<ZZ_pEX> &auxValuesT,
                         const Vec<ZZ_pEX> &values);

    void computAuxValues(Vec<ZZ_pE> &auxValuesS,
                         Vec<ZZ_pE> &auxValuesT,
                         const Vec<ZZ_pE> &values);

    long modulusDegree;
    ZZ_pE a;
    ZZ_pE b;
    ZZ_pEX psi2SquareEX;
    ZZ_pE psi2SquareE;
    ZZ_pEX psi2;
};

#endif /* DIVISIONPOLY_H */

