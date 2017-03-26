

#ifndef ISOGENYGRAPH_H
#define ISOGENYGRAPH_H

#include <NTL/ZZ_pEXFactoring.h>

using namespace NTL;

class IsogenyGraph {
public:

    IsogenyGraph();
    virtual ~IsogenyGraph();

    void getRandomSuperSingular(ZZ_pE &a1,
                                ZZ_pE &b1,
                                const ZZ_pE &a,
                                const ZZ_pE &b, long m);
    ZZ_pE getJInvariant(const ZZ_pE &a, const ZZ_pE &b);
    ZZ_pEX getModularPolynomial(const ZZ_pE &j);
    bool isSupersingular(const ZZ_pE &a, const ZZ_pE &b);

private:

    ZZ_p a;
    ZZ_p b;
    ZZ_p c;
    ZZ_p d;
    ZZ_p e;
};

#endif /* ISOGENYGRAPH_H */

