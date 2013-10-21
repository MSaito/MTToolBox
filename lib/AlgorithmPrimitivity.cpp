#include <sstream>
#include <NTL/GF2XFactoring.h>
#include <NTL/ZZ.h>
#include <NTL/vector.h>
#include <MTToolBox/AlgorithmPrimitivity.hpp>
#include <MTToolBox/period.hpp>


namespace MTToolBox {
    using namespace std;
    using namespace NTL;

    AlgorithmPrimitivity::AlgorithmPrimitivity(const char * prime_factors[])
    {
        mersenne = false;
        stringstream ss;
        long length = 0;
        for (int i = 0; prime_factors[i]; i++) {
            length = i + 1;
        }
        this->primes = new Vec<ZZ>;
        this->primes->SetLength(length);
        for (int i = 0; prime_factors[i]; i++) {
            ss << prime_factors[i];
            ss << " ";
            ZZ w;
            ss >> w;
            (*primes)[i] = w;
        }
    }

    bool AlgorithmPrimitivity::operator()(int max_degree,
                                          const NTL::GF2X& poly) const
    {
        if (max_degree != deg(poly)) {
            return false;
        }
        if (!IterIrredTest(poly)) {
            return false;
        }
        if (mersenne) {
            return true;
        }
        if (isPrime(poly, max_degree, *primes)) {
            return true;
        } else {
            return false;
        }
    }

    const AlgorithmPrimitivity MersennePrimitivity;
}
