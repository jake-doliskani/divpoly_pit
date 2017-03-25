#include "IsogenyGraph.h"

IsogenyGraph::IsogenyGraph() {
}

IsogenyGraph::~IsogenyGraph() {
}

ZZ_pE IsogenyGraph::getJInvariant(const ZZ_pE &a, const ZZ_pE &b) {
	return 1728 * (4 * a * a * a) / (4 * a * a * a + 27 * b * b);
}

ZZ_pEX IsogenyGraph::getModularPolynomial(const ZZ_pE &j) {
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

void IsogenyGraph::getRandomSuperSingular(ZZ_pE &a1,
										  ZZ_pE &b1,
										  const ZZ_pE &a,
										  const ZZ_pE &b, long m) {
	ZZ_pE jInv = getJInvariant(a, b);
	ZZ_pEX phi2 = getModularPolynomial(jInv);
	Vec<ZZ_pE> rootVec;

	ZZ_pEX tempX;
	SetX(tempX);

	for (long i = 0; i <= m; i++) {
		rootVec = FindRoots(phi2);
		phi2 = getModularPolynomial(rootVec[0]) / (tempX - jInv);
		jInv = rootVec[0];
	}

	a1 = 3 * jInv / (1728 - jInv);
	b1 = 2 * jInv / (1728 - jInv);
}

bool IsogenyGraph::isSupersingular(const ZZ_pE &a, const ZZ_pE &b) {

	ZZ_pE j = getJInvariant(a, b);
	ZZ_pEX phi2 = getModularPolynomial(j);

	Vec<pair_ZZ_pEX_long> temp = CanZass(phi2);
	if (temp[0].b > 1 || temp[1].b > 1)
		return false;

	long m = NumBits(ZZ_p::modulus()) + 1;

	Vec<ZZ_pE> jInvVec;
	jInvVec.SetLength(3);
	for (int i = 0; i < 3; i++)
		jInvVec[i] = j;
	
	Vec<ZZ_pE> rootVec;
	rootVec.SetLength(3);
	for (int i = 0; i < 3; i++)
		rootVec[i] = -coeff(temp[i].a, 0);
	
	ZZ_pEX tempX;
	SetX(tempX);

	Vec<ZZ_pE> tempRoots;
	for (long i = 0; i <= m; i++) {
		for (long j = 0; j < 3; j++) {
			phi2 = getModularPolynomial(rootVec[j]) / (tempX - jInvVec[j]);
			if (ProbIrredTest(phi2))
				return false;
			
			tempRoots = FindRoots(phi2);
			jInvVec[j] = rootVec[j];
			rootVec[j] = tempRoots[0];
		}
	}

	return true;
}
