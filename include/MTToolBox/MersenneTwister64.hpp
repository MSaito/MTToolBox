#ifndef MTTOOLBOX_MERSENNETWISTER64_HPP
#define MTTOOLBOX_MERSENNETWISTER64_HPP

/**
 * @file MersenneTwister64.hpp
 *
 *\japanese
 * @brief 64bit MersenneTwister generator
 *\endjapanese
 *
 *\english
 * @brief 64 bit MersenneTwister generator
 *\endenglish
 *
 *
 * Original
 * @author Makoto Matsumoto (Hiroshima University)
 * @author Takuji Nishimura (Yamagata University)
 * Copyright (c) 2004 Makoto Matsumoto and Takuji Nishimura.
 * All rights reserved.
 *
 * Modified
 * @author Mutsuo Saito (Manieth Corp.)
 * @author Makoto Matsumoto (Hiroshima University)
 * Copyright (c) 2015 Mutsuo Saito, Makoto Matsumoto and Hiroshima
 * University. All rights reserved.
 *
 * The 3-clause BSD License is applied to this software, see
 * LICENSE.txt
 */

#include <stdint.h>
#include <inttypes.h>
#include <string>
#include <MTToolBox/ParameterGenerator.hpp>
#if defined(DEBUG)
#include <iostream>
#include <iomanip>
#endif

namespace MTToolBox {
    /**
     *\japanese
     * 64 bit Mersenne Twister 疑似乱数生成器
     *\endjapanese
     *
     *\english
     * Mersenne Twister pseudo random number generator.
     *
     *
     *\endenglish
     */
    class MersenneTwister64 : public ParameterGenerator {
    public:
        /**
         *\japanese
         * コンストラクタ
         *\endjapanese
         *
         *\english
         * Constructor without arguments.
         *\endenglish
         */
        MersenneTwister64() {
            mt = new uint64_t[N];
            seed(UINT64_C(19650218));
        }

        /**
         *\japanese
         * コンストラクタ
         * @param[in] value 初期化の種
         *\endjapanese
         *
         *\english
         * Constructor with a integer seed
         * @param[in] value Seed of initialization
         *\endenglish
         */
        MersenneTwister64(uint64_t value) {
            mt = new uint64_t[N];
            seed(value);
        }

        /**
         *\japanese
         * コンストラクタ
         * @param[in] value 初期化の種（文字列）
         *\endjapanese
         *
         *\english
         * Constructor with a string seed
         * @param[in] value Seed of initialization
         *\endenglish
         */
        MersenneTwister64(const std::string& value) {
            mt = new uint64_t[N];
            seed(value);
        }

        /**
         *\japanese
         * コンストラクタ
         * @param[in] value 初期化の種（配列）
         * @param[in] size 配列の長さ
         *\endjapanese
         *
         *\english
         * Constructor with an array of integer
         * @param[in] value Seed of initialization
         * @param[in] size length of \b value
         *\endenglish
         */
        MersenneTwister64(const uint64_t *value, int size) {
            mt = new uint64_t[N];
            seed(value, size);
        }

        /**
         *\japanese
         * デストラクタ
         *\endjapanese
         *
         *\english
         * Destructor
         *\endenglish
         */
        ~MersenneTwister64() {
            delete[] mt;
        }

        /**
         *\japanese
         * 64bit整数による初期化
         * @param[in] value 初期化の種
         *\endjapanese
         *
         *\english
         * Initialization by 64-bit integer
         * @param[in] value Seed of initialization
         *\endenglish
         */
        void seed(uint64_t value) {
            mt[0] = value;
            for (mti = 1; mti < N; mti++) {
                mt[mti] = mti
                    + UINT64_C(6364136223846793005)
                    * (mt[mti-1] ^ (mt[mti-1] >> 62));
            }
            mti = 0;
        }

        /**
         *\japanese
         * string文字列による初期化
         * @param[in] value 初期化の種
         *\endjapanese
         *
         *\english
         * Initialization by a string
         * @param[in] value Seed of initialization
         *\endenglish
         */
        void seed(const std::string& value) {
            seed_array<char>(value.c_str(), static_cast<int>(value.size()));
        }

        /**
         *\japanese
         * 符号なし整数配列による初期化
         * @param[in] value 初期化の種
         * @param[in] key_length 配列の長さ
         *\endjapanese
         *
         *\english
         * Initialization by an array of unsigned integers.
         * @param[in] value Seed of initialization
         * @param[in] key_length Length of \b value
         *\endenglish
         */
        void seed(const uint64_t *value, int key_length) {
            seed_array<uint64_t>(value, key_length);
        }

        /**
         *\japanese
         * T 型配列による初期化
         * @tparam T 配列の要素の型
         * @param[in] init_key 初期化の種
         * @param[in] key_length 文字列の長さ
         *\endjapanese
         *
         *\english
         * Initialization by an array of T type
         * @tparam T type of element of \b value
         * @param[in] init_key Seed of initialization
         * @param key_length Length of \b value
         *\endenglish
         */
        template<class T> void seed_array(const T *init_key,
                                          int key_length) {
            int i, j;
            uint64_t k;
            seed(UINT64_C(19650218));
            i = 1;
            j = 0;
            k = (N > key_length ? N : key_length);
            for (; k; k--) {
                mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 62))
                                  * UINT64_C(3935559000370003845)))
                    + init_key[j] + (uint64_t)j; /* non linear */
                i++;
                j++;
                if (i>=N) {
                    mt[0] = mt[N-1];
                    i = 1;
                }
                if (j >= key_length) {
                    j=0;
                }
            }
            for (k=N-1; k; k--) {
                mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 62))
                                  * UINT64_C(2862933555777941757)))
                    - (uint64_t)i; /* non linear */
                i++;
                if (i>=N) {
                    mt[0] = mt[N-1];
                    i=1;
                }
            }
            /* MSB is 1; assuring non-zero initial array */
            mt[0] = UINT64_C(1) << 63;
            mti = 0;
        }

        /**
         *\japanese
         * 疑似乱数を生成する
         * @return １個の32bit符号なし整数
         *\endjapanese
         *
         *\english
         * Generates pseudo random number
         * @return a 32-bit unsigned integer
         *\endenglish
         */
        uint32_t getUint32() {
            return next() >> 32;
        }

        /**
         *\japanese
         * 疑似乱数を生成する
         * @return １個の64bit符号なし整数
         *\endjapanese
         *
         *\english
         * Generates pseudo random number
         * @return a 64-bit unsigned integer
         *\endenglish
         */
        uint64_t getUint64() {
            return next();
        }

        /**
         *\japanese
         * 疑似乱数を生成する
         * @return １個の64bit符号なし整数
         *\endjapanese
         *
         *\english
         * Generates pseudo random number
         * @return a 64-bit unsigned integer
         *\endenglish
         */
        uint64_t next() {
            const uint64_t UPPER_MASK = UINT64_C(0xFFFFFFFF80000000);
            const uint64_t LOWER_MASK = UINT64_C(0x000000007FFFFFFF);
            const uint64_t mag01[2] = {0x0UL, UINT64_C(0xB5026F5AA96619E9)};
            uint64_t x;

            x = (mt[mti] & UPPER_MASK) | (mt[(mti + 1) % N] & LOWER_MASK);
            mt[mti] = mt[(mti + M) % N] ^ (x >> 1)
                ^ mag01[(int)(x & UINT64_C(1))];
            x = mt[mti];
            mti = (mti + 1) % N;
            return temper(x);
        }

        /**
         *\japanese
         * 状態空間のビットサイズである 19937 を返す。
         * @return 常に 19937 を返す
         *\endjapanese
         *
         *\english
         * Returns 19937, which is size of internal state.
         * @return Always 19937
         *\endenglish
         */
        int bitSize() const {
            return 19937;
        }
    private:
        enum {N = 312, M = 156};
        uint64_t *mt;    /* the array for the state vector  */
        int mti;
        uint64_t temper(uint64_t y) {
            y ^= (y >> 29) & UINT64_C(0x5555555555555555);
            y ^= (y << 17) & UINT64_C(0x71D67FFFEDA60000);
            y ^= (y << 37) & UINT64_C(0xFFF7EEE000000000);
            y ^= (y >> 43);
            return y;
        }
    };
}
#endif // MTTOOLBOX_MERSENNETWISTER64_HPP
