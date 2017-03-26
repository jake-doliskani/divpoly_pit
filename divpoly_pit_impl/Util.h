#ifndef UTIL_H
#define UTIL_H

#include <NTL/ZZ_pX.h>

using namespace NTL;

class Util {
public:

    Util();
    virtual ~Util();

    /**
     * @return the system time in milliseconds
     */
    long getTimeMillis();

    long findOrder(long a, long modulus);

private:

    class Factor {
    public:
        long p;
        long e;
    };

    void factorNaive(Vec<Factor> &factors, long a);

    long EulerTotient(long a);
};

#endif /* UTIL_H */

