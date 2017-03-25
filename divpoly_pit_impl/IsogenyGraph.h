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

};

#endif /* ISOGENYGRAPH_H */

