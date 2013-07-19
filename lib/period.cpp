#include <NTL/GF2XFactoring.h>
#include <NTL/ZZ.h>
#include <NTL/vector.h>
#include <MTToolBox/period.hpp>


namespace MTToolBox {
    using namespace std;
    using namespace NTL;

    static const uint32_t mersenne_exponent[] =
    {2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279,
     2203, 2281, 3217, 4253, 4423, 9689, 9941, 11213, 19937, 21701,
     23209, 44497, 86243, 110503, 132049, 216091, 756839, 859433,
     1257787, 1398269, 2976221, 3021377, 6972593, 13466917,
     20996011, 25964951, -1};

    void
    minpoly(NTL::GF2X& poly, U128Generator& generator, int pos)
    {
        Vec<GF2> v;
        int size = generator.bitSize();
        v.SetLength(2 * size);
        int index;
        if (pos < 64) {
            index = 0;
        } else {
            index = 1;
            pos = pos - 64;
        }
        for (int i = 0; i < 2 * size; i++) {
            uint128_t u = generator.generate();
            v[i] = (u.u64[index] >> pos) & 1;
        }
        MinPolySeq(poly, v, size);
    }

    bool
    isIrreducible(NTL::GF2X& poly)
    {
        return static_cast<bool>(IterIrredTest(poly));
    }

    bool
    isPrime(NTL::GF2X& poly)
    {
        if (!isIrreducible(poly)) {
            return false;
        }
        long degree = deg(poly);
        for (int i = 0; mersenne_exponent[i] > 0; i++) {
            if (degree == mersenne_exponent[i]) {
                return true;
            }
        }
        return false;
    }

    bool
    isPrime(NTL::GF2X& poly, NTL::Vec<ZZ>& prime_factors)
    {
        if (!isIrreducible(poly)) {
            return false;
        }
        long degree = deg(poly);
        long len = primes.length();
        ZZ period(2);
        period = Power(perid, degree);
        period -= 1;
        for (long i = 0; i < len; i++) {
            ZZ p = prime_factors[i];
            ZZ pow = period / p;
        }
        return false;
    }
}
