#include <MTToolBox/period.hpp>
#include <MTToolBox/AlgorithmPrimitivity.hpp>
#include <UnitTest++.h>
#include <NTL/ZZ.h>

using namespace MTToolBox;
using namespace NTL;
using namespace std;

void ctozz(ZZ& w, const char * large_int) {
    stringstream ss;
    ss << large_int << endl;
    ss >> w;
}

SUITE(PRIME_FACTORS) {
    TEST(factors)
    {
        const char ** primes[] = {
            prime_factors2_128_1,
            prime_factors2_160_1,
            prime_factors2_192_1,
            prime_factors2_224_1,
            prime_factors2_256_1,
            prime_factors2_288_1,
            prime_factors2_320_1,
            prime_factors2_352_1,
            prime_factors2_384_1,
            prime_factors2_416_1,
            prime_factors2_448_1,
            prime_factors2_480_1,
            prime_factors2_512_1,
            prime_factors2_544_1,
            NULL};
        ZZ two32;
        two32 = power2_ZZ(32);
        for (int i = 4; i <= 17; i++) {
            ZZ w = two32;
            w = power(w, i);
            w -= 1;
            const char ** p = primes[i - 4];
            for (int j = 0; p[j] != NULL; j++) {
                ZZ div;
                ctozz(div, p[j]);
                ZZ r;
                r = w % div;
                CHECK(IsZero(r));
            }
        }
    }
}
