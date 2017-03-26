#ifndef CYCLOMOD_H
#define CYCLOMOD_H

#include <NTL/ZZ_pEX.h>

using namespace NTL;

class CycloMod {
public:

    CycloMod(long modulusDegree);
    virtual ~CycloMod();

    void sqrMod(ZZ_pEX &x, const ZZ_pEX& f);

    void reduce(ZZ_pEX &x, const ZZ_pEX &a);

    void mulMod(ZZ_pEX &x, const ZZ_pEX &f, const ZZ_pEX &g);

    void invMod(ZZ_pEX &x, const ZZ_pEX &f);

    //    void mul(ZZ_pE &x, const ZZ_pE &a, const ZZ_pE &b);

private:

    long modulusDegree;
};

#endif /* CYCLOMOD_H */

