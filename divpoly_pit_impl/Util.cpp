
#include "Util.h"
#include <chrono>
#include <NTL/FacVec.h>

using namespace std;

Util::Util() {
}

Util::~Util() {
}

long Util::getTimeMillis() {
    auto now = chrono::system_clock::now();
    auto now_ms = chrono::time_point_cast<chrono::milliseconds>(now);
    auto value = now_ms.time_since_epoch();
    long duration = value.count();
    return duration;
}

long Util::findOrder(long a, long modulus) {
	if (GCD(a, modulus) > 1) {
		return -1;
	}
	
	long order = EulerTotient(modulus);
    Vec<Factor> factors;
	factorNaive(factors, order);
	
	long temp = 1;
	for (long i = 0; i < factors.length(); i++) {
		while (temp == 1 && factors[i].e > 0) {
			order /= factors[i].p;
			factors[i].e--;
			temp = PowerMod(a, order, modulus);
		}

		if (temp != 1) {
			order *= factors[i].p;
			temp = 1;
		}
	}

	return order;
}

void Util::factorNaive(Vec<Factor>& factors, long a) {
	long p = 2;
	
	while (a > 1) {
		long i = 0;
		
		while (a % p == 0) {
			a = a / p;
			i++;
		}
		
		if (i > 0) {
			Factor factor;
			factor.p = p;
			factor.e = i;
			factors.append(factor);
		}
		
		p = NextPrime(p + 1);
	}
}

long Util::EulerTotient(long a) {
    Vec<Factor> factors;
	factorNaive(factors, a);
	
	long result = 1;
	for (long i = 0; i < factors.length(); i++) {
		result *= power_long(factors[i].p, factors[i].e - 1);
		result *= factors[i].p - 1;
	}
	
	return result;
}