import time


def readPrimes(filename):
	f = open(filename, "r")
	primes = []
	for i in range(0, 10):
		f.readline()
		primes.append(int(f.readline()))
		f.readline()
		f.readline()

	f.close()
	return primes


def buildExtensionAndCurve(p):
	Fpp.<i> = GF(p^2, modulus = x^2 + 1)
	E = EllipticCurve(Fpp, [Fpp.random_element(), Fpp.random_element()])
	return Fpp, E


def currentTimeMillis():
	return int(round(time.time() * 1000))


def main():
	primes = readPrimes("supersingular_params")
	f = open("sage_times", "w", 0)
	f.write("Numbits  time\n")
	for i in range(0, len(primes)):
		Fpp, E = buildExtensionAndCurve(primes[i])
		time = currentTimeMillis()
		E._pari_().ellsea()
		time = currentTimeMillis() - time
		f.write(str(ZZ(primes[i]).nbits()) + "     " + str(time) + "\n")

	f.close()
