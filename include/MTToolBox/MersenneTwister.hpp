#ifndef MTTOOLBOX_MERSENNETWISTER_HPP
#define MTTOOLBOX_MERSENNETWISTER_HPP

/**
 * @file MersenneTwister.hpp
 *
 * @brief MersenneTwister generator
 *
 * @author Mutsuo Saito (Hiroshima University)
 * @author Makoto Matsumoto (Hiroshima University)
 *
 * Copyright (c) 2013 Mutsuo Saito, Makoto Matsumoto and Hiroshima
 * University. All rights reserved.
 *
 * The 3-clause BSD License is applied to this software, see
 * LICENSE.txt
 */

#include <stdint.h>
#include <inttypes.h>
#include <string>
#include <MTToolBox/AbstractGenerator.hpp>

namespace MTToolBox {
    /**
     * Mersenne Twister 疑似乱数生成器
     *
     * NOTE: この実装では、AbstractGenerator の仕様にあわせて、
     * 一度の呼び出しで１個の疑似乱数を生成するようにしてある。
     * そうするとまとめて作る場合よりも生成速度が落ちるので、
     * 気休め程度でも速度を向上させるために、状態空間より大きな配列を使用して
     * 剰余演算をビットマスクで済ませるようにしてある。
     */
    class MersenneTwister : public AbstractGenerator<uint32_t> {
    public:
        /**
         * コンストラクタ
         */
        MersenneTwister() {
            mt = new uint32_t[LARGE_N];
            seed(5489);
        }

        /**
         * コンストラクタ
         * @param[in] value 初期化の種
         */
        MersenneTwister(uint32_t value) {
            mt = new uint32_t[LARGE_N];
            seed(value);
        }

        /**
         * コンストラクタ
         * @param[in] value 初期化の種（文字列）
         */
        MersenneTwister(const std::string& value) {
            mt = new uint32_t[LARGE_N];
            seed(value);
        }

        /**
         * コンストラクタ
         * @param[in] value 初期化の種（配列）
         * @param[in] size 配列の長さ
         */
        MersenneTwister(const uint32_t *value, int size) {
            mt = new uint32_t[LARGE_N];
            seed(value, size);
        }

        /**
         * デストラクタ
         */
        ~MersenneTwister() {
            delete[] mt;
        }

        /**
         * 32bit整数による初期化
         * @param[in] value 初期化の種
         */
        void seed(uint32_t value) {
            mt[0] = value;
            for (mti = 1; mti < N; mti++) {
                mt[mti] =
                    (UINT32_C(1812433253)
                     * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + mti);
            }
        }

        /**
         * string文字列による初期化
         * @param[in] value 初期化の種
         */
        void seed(const std::string& value) {
            seed_array<char>(value.c_str(), static_cast<int>(value.size()));
        }

        /**
         * char 文字列による初期化
         * @param[in] value 初期化の種
         * @param[in] key_length 文字列の長さ
         */
        void seed(const uint32_t *value, int key_length) {
            seed_array<uint32_t>(value, key_length);
        }

        /**
         * T 型配列による初期化
         * @tparam T 配列の要素の型
         * @param[in] value 初期化の種
         * @param[in] key_length 文字列の長さ
         */
        template<class T> void seed_array(const T *value, int key_length) {
            int i, j, k;
            seed(UINT32_C(19650218));
            i = 1;
            j = 0;
            if (N > key_length) {
                k = N;
            } else {
                k = key_length;
            }
            for (; k > 0; k--) {
                mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30))
                                  * UINT32_C(1664525)))
                    + static_cast<uint32_t>(value[j]) + j; /* non linear */
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
                mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30))
                                  * UINT32_C(1566083941)))
                    - i; /* non linear */
                i++;
                if (i>=N) {
                    mt[0] = mt[N-1];
                    i=1;
                }
            }
            /* MSB is 1; assuring non-zero initial array */
            mt[0] = UINT32_C(0x80000000);
            mti = N;
        }

        /**
         * 疑似乱数を生成する
         * @return １個の32bit符号なし整数
         */
        uint32_t generate() {
            return next();
        }

        /**
         * 疑似乱数を生成する
         * @return １個の32bit符号なし整数
         */
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

        /**
         * 状態空間のビットサイズである 19937 を返す。
         *
         * 配列は状態空間より大きく取ってあるが、それは速度向上のために
         * 冗長になっているだけで、状態空間の大きさはあくまでも19937である。
         * @return 常に 19937 を返す
         */
        int bitSize() const {
            return 19937;
        }
    private:
        enum {LARGE_N = 1024, N = 624, M = 397};
        uint32_t *mt;    /* the array for the state vector  */
        int mti;
        uint32_t temper(uint32_t y) {
            y ^= (y >> 11);
            y ^= (y << 7) & UINT32_C(0x9d2c5680);
            y ^= (y << 15) & UINT32_C(0xefc60000);
            y ^= (y >> 18);
            return y;
        }
    };
}
#endif // MTTOOLBOX_MERSENNETWISTER_HPP
