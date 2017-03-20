#ifndef CYCLOMOD_H
#define CYCLOMOD_H

#include <NTL/ZZ_pEX.h>

using namespace NTL;


class CycloMod {
    
public:
    
    CycloMod(long modulusDegree);
    virtual ~CycloMod();
    
    void sqrMod(ZZ_pEX& result, const ZZ_pEX& f);
    
    void reduce(ZZ_pEX &result, const ZZ_pEX &a);

    void mulMod(ZZ_pEX &result, const ZZ_pEX &f, const ZZ_pEX &g);

    void powerMod(ZZ_pEX &result, const ZZ_pEX &f, const ZZ &e);
    
private:

    long optimalWinSize(long n);
    
    long modulusDegree;
    
};

#endif /* CYCLOMOD_H */

