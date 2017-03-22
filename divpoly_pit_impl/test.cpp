#include <cstdlib>
#include "DivisionPoly.h"
#include "CycloMod.h"
#include "Util.h"
#include "DivpolyPit.h"
#include <iostream>
#include <NTL/ZZ_pXFactoring.h>
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

void testSupersingular() {
    vector<TestCase> testCases;
    char fileName[100] = "supersingular_params";
    readTestCases(testCases, fileName);
    
    cout << "numBits time type" << endl;
    for (TestCase testCase : testCases) {
        testPIT(testCase);
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
    testOrdinary();
    
    return 0;
}

