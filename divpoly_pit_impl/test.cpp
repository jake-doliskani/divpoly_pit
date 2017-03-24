#include <cstdlib>
#include "DivisionPoly.h"
#include "CycloMod.h"
#include "Util.h"
#include "DivpolyPit.h"
#include <iostream>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pEXFactoring.h>
#include <NTL/ZZX.h>
#include <vector>

using namespace std;

class TestCase {
public:
    long numBits;
    ZZ p;
    ZZ a0;
    ZZ a1;
    ZZ b0;
    ZZ b1;
};

void readTestCases(vector<TestCase> &testCases, char* fileName) {
    FILE *input = fopen(fileName, "r");

    while (true) {
        char str1[1000];
        char str2[1000];
        TestCase testCase;

        if (fscanf(input, "%li\n", &testCase.numBits) == EOF)
            break;

        if (fscanf(input, "%s\n", str1) == EOF)
            break;
        testCase.p = to_ZZ(str1);

        if (fscanf(input, "%s%*[*i + ]%s\n", str1, str2) == EOF)
            break;
        testCase.a1 = to_ZZ(str1);
        testCase.a0 = to_ZZ(str2);

        if (fscanf(input, "%s%*[*i + ]%s\n", str1, str2) == EOF)
            break;
        testCase.b1 = to_ZZ(str1);
        testCase.b0 = to_ZZ(str2);

        testCases.push_back(testCase);
    }

    fclose(input);
}

void testPIT(TestCase testCase) {
    ZZ_p::init(testCase.p);

    ZZ_pX f;
    clear(f);
    SetCoeff(f, 2, 1);
    SetCoeff(f, 0, 1);
    ZZ_pE::init(f);

    ZZ_pE a, b;
    SetCoeff(a._ZZ_pE__rep, 0, to_ZZ_p(testCase.a0));
    SetCoeff(a._ZZ_pE__rep, 1, to_ZZ_p(testCase.a1));
    SetCoeff(b._ZZ_pE__rep, 0, to_ZZ_p(testCase.b0));
    SetCoeff(b._ZZ_pE__rep, 1, to_ZZ_p(testCase.b1));

    DivpolyPit divpolyPit;
    Util util;
    long time = util.getTimeMillis();
    bool isSupersingular = divpolyPit.isSuperSingular(a, b);
    time = util.getTimeMillis() - time;

    cout << testCase.numBits << " " << time << " ";
    if (isSupersingular)
        cout << "supersingular" << endl;
    else
        cout << "ordinary" << endl;
}

//////////////////////////////////////////////////

ZZ_pE getJInvariant(const ZZ_pE &a, const ZZ_pE &b) {
    return 1728 * (4 * a * a * a) / (4 * a * a * a + 27 * b * b);
}

ZZ_pEX getModularPolynomial(const ZZ_pE &j) {
    ZZ_p a = to_ZZ_p(to_ZZ("1488"));
    ZZ_p b = to_ZZ_p(to_ZZ("162000"));
    ZZ_p c = to_ZZ_p(to_ZZ("40773375"));
    ZZ_p d = to_ZZ_p(to_ZZ("8748000000"));
    ZZ_p e = to_ZZ_p(to_ZZ("157464000000000"));

    ZZ_pE a2 = -j * j + a * j - b;
    ZZ_pE a1 = a * j * j + c * j + d;
    ZZ_pE a0 = j * j * j - b * j * j + d * j - e;

    ZZ_pEX result;
    clear(result);
    SetCoeff(result, 3, 1);
    SetCoeff(result, 2, a2);
    SetCoeff(result, 1, a1);
    SetCoeff(result, 0, a0);

    return result;
}

bool testIsogeyGraph(TestCase testCase) {
    ZZ_p::init(testCase.p);

    ZZ_pX f;
    clear(f);
    SetCoeff(f, 2, 1);
    SetCoeff(f, 0, 1);
    ZZ_pE::init(f);

    ZZ_pE a, b;
    SetCoeff(a._ZZ_pE__rep, 0, to_ZZ_p(testCase.a0));
    SetCoeff(a._ZZ_pE__rep, 1, to_ZZ_p(testCase.a1));
    SetCoeff(b._ZZ_pE__rep, 0, to_ZZ_p(testCase.b0));
    SetCoeff(b._ZZ_pE__rep, 1, to_ZZ_p(testCase.b1));

    ZZ_pE j = getJInvariant(a, b);
    ZZ_pEX phi2 = getModularPolynomial(j);

    Vec<pair_ZZ_pEX_long> tempJVec = CanZass(phi2);
    if (tempJVec.length() < 3)
        return false;

    long m = testCase.numBits + 1;

    Vec<ZZ_pEX> phi2Vec;
    phi2Vec.SetLength(3);

    Vec<ZZ_pE> jVec;
    jVec.SetLength(3);
    for (long i = 0; i < 3; i++)
        jVec[i] = -coeff(tempJVec[i].a, 0);

    ZZ_pEX tempX;
    SetX(tempX);

    for (long i = 0; i <= m; i++) {

        for (long j = 0; j < 3; j++) {
            phi2Vec[j] = getModularPolynomial(jVec[j]) / (tempX - jVec[j]);
        }

        tempJVec.SetLength(0);
        for (long j = 0; j < 3; j++) {
            tempJVec = CanZass(phi2Vec[j]);
            if (tempJVec.length() < 2)
                return false;

            jVec[j] = -coeff(tempJVec[0].a, 0);
        }
    }

    return true;
}





void testSupersingular() {
    vector<TestCase> testCases;
    char fileName[100] = "supersingular_params";
    readTestCases(testCases, fileName);

    cout << "numBits time type" << endl;
    for (TestCase testCase : testCases) {
        Util util;
        long start = util.getTimeMillis();
        testIsogeyGraph(testCase);
        cout << testCase.numBits << " " << util.getTimeMillis() - start << endl;
        //        testPIT(testCase);
    }
}

void testOrdinary() {
    vector<TestCase> testCases;
    char fileName[100] = "ordinary_params";
    readTestCases(testCases, fileName);

    cout << "numBits time type" << endl;
    for (TestCase testCase : testCases) {
        testPIT(testCase);
    }
}

int main(int argc, char** argv) {

    testSupersingular();
//    testOrdinary();

    return 0;
}

