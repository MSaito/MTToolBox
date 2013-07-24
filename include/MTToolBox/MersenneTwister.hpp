#ifndef MT19937_HPP
#define MT19937_HPP

/* This macro should be defined compiler option or .c .cpp file */
#ifndef __STDC_CONSTANT_MACROS
#define __STDC_CONSTANT_MACROS 1
#endif

/**
 * @file all_in_one.hpp
 *
 * @brief MersenneTwister generator for generating random parameters.
 *
 * @author Mutsuo Saito (Hiroshima University)
 * @author Makoto Matsumoto (Hiroshima University)
 *
 * Copyright (c) 2010 Mutsuo Saito, Makoto Matsumoto and Hiroshima
 * University. All rights reserved.
 * Copyright (c) 2011 Mutsuo Saito, Makoto Matsumoto, Hiroshima
 * University and University of Tokyo. All rights reserved.
 *
 * The new BSD License is applied to this software, see LICENSE.txt
 */

#include <stdint.h>
#include <inttypes.h>
#include <string>

namespace MTToolBox {
    class MersenneTwister {
    public:
        MersenneTwister() {
            mt = new uint32_t[LARGE_N];
            reseed(5489);
        }
        MersenneTwister(uint32_t seed) {
            mt = new uint32_t[LARGE_N];
            reseed(seed);
        }
        MersenneTwister(const std::string& seed) {
            mt = new uint32_t[LARGE_N];
            reseed(seed);
        }
        MersenneTwister(const uint32_t *seed, int size) {
            mt = new uint32_t[LARGE_N];
            reseed(seed, size);
        }
        ~MersenneTwister() {
            delete[] mt;
        }
        void reseed(uint32_t seed) {
            mt[0] = seed;
            for (mti = 1; mti < N; mti++) {
                mt[mti] =
                    (1812433253UL * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + mti);
            }
        }
        void reseed(const std::string& seed) {
            reseed<char>(seed.c_str(), static_cast<int>(seed.size()));
        }
        void reseed(const uint32_t *seed, int key_length) {
            reseed<uint32_t>(seed, key_length);
        }
        template<class T> void reseed(const T *seed, int key_length) {
            int i, j, k;
            reseed(19650218UL);
            i = 1;
            j = 0;
            if (N > key_length) {
                k = N;
            } else {
                k = key_length;
            }
            for (; k > 0; k--) {
                mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
                    + static_cast<uint32_t>(seed[j]) + j; /* non linear */
                i++;
                j++;
                if (i>=N) {
                    mt[0] = mt[N-1];
                    i=1;
                }
                if (j>=key_length) {
                    j=0;
                }
            }
            for (k=N-1; k; k--) {
                mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
                    - i; /* non linear */
                i++;
                if (i>=N) {
                    mt[0] = mt[N-1];
                    i=1;
                }
            }
            /* MSB is 1; assuring non-zero initial array */
            mt[0] = 0x80000000UL;
            mti = N;
        }
        uint32_t next() {
            using namespace std;

            const uint32_t LARGE_MASK = UINT32_C(0x3ff);
            const uint32_t UPPER_MASK = UINT32_C(0x80000000);
            const uint32_t LOWER_MASK = UINT32_C(0x7fffffff);
            const uint32_t mag01[2] = {0x0UL, UINT32_C(0x9908b0df)};
            uint32_t y;

            y = (mt[(LARGE_N - N + mti) & LARGE_MASK] & UPPER_MASK)
                | (mt[(LARGE_N - N + mti + 1) & LARGE_MASK] & LOWER_MASK);
            mt[mti] = mt[(LARGE_N - N + M + mti) & LARGE_MASK]
                ^ (y >> 1) ^ mag01[y & 0x1UL];
            y = temper(mt[mti]);
            mti = (mti + 1) & LARGE_MASK;
            return y;
        }
    private:
        enum {LARGE_N = 1024, N = 624, M = 397};
        uint32_t *mt;    /* the array for the state vector  */
        int mti;
        uint32_t temper(uint32_t y) {
            y ^= (y >> 11);
            y ^= (y << 7) & 0x9d2c5680UL;
            y ^= (y << 15) & 0xefc60000UL;
            y ^= (y >> 18);
            return y;
        }
    };
    extern MersenneTwister MT;
}
#endif
