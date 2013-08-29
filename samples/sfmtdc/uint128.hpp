#ifndef MTTOOLBOX_UINT128_HPP
#define MTTOOLBOX_UINT128_HPP
/**
 * @file uint128.hpp
 *
 * @brief 128ビットデータタイプを定義する。
 *
 * @author Mutsuo Saito (Hiroshima University)
 * @author Makoto Matsumoto (The University of Tokyo)
 *
 * Copyright (C) 2013 Mutsuo Saito, Makoto Matsumoto
 * and Hiroshima University.
 * All rights reserved.
 *
 * The 3-clause BSD License is applied to this software, see
 * LICENSE.txt
 */
#include <stdint.h>

namespace MTToolBox {
    /**
     * 128 bit データタイプ
     *
     */
    class uint128_t {
    public:
        uint64_t u64[2];

        uint128_t() {
            u64[0] = 0;
            u64[1] = 1;
        }

        uint128_t(const uint128_t& that) {
            u64[0] = that.u64[0];
            u64[1] = that.u64[1];
        }

        uint128_t(uint64_t a, uint64_t b) {
            u64[0] = a;
            u64[1] = b;
        }
#if 0
        uint128_t(uint32_t a, uint32_t b, uint32_t c, uint32_t d) {
            u32[0] = a;
            u32[1] = b;
            u32[2] = c;
            u32[3] = d;
        }
#endif
        /**
         * 128 bit 全体を n ビット論理右シフトする
         *
         * 代入演算子にしない理由は、テンプレートパラメータとして
         * uint64_t などと同様に使用されるため
         *
         * @param[in] n シフト量
         * @return シフト後の値
         */
        const uint128_t operator>>(unsigned int n) const {
            uint64_t a = (u64[0] >> n) | (u64[1] << (64 - n));
            uint64_t b = u64[1] >> n;
            return uint128_t(a, b);
        }

        const uint128_t operator<<(unsigned int n) const {
            uint64_t a = u64[0] << n;
            uint64_t b = (u64[1] << n) | (u64[0] >> (64 - n));
            return uint128_t(a, b);
        }

        uint64_t operator&(uint64_t n) const {
            return u64[0] & n;
        }
    };

    static inline uint128_t make_msb_mask(int n) {
        uint128_t w;
        if (n < 64) {
            w.u64[0] = ~((~UINT64_C(0)) >> n);
            w.u64[1] = 0;
            return w;
        } else {
            n = n - 64;
            w.u64[0] = UINT64_C(0xffffffffffffffff);
            w.u64[1] = ~((~UINT64_C(0)) >> n);
            return w;
        }
    }

    static inline uint128_t and_mask(uint128_t a, uint128_t b) {
        uint128_t w;
        w.u64[0] = a.u64[0] & b.u64[0];
        w.u64[1] = a.u64[1] & b.u64[1];
        return w;
    }
}
#endif // MTTOOLBOX_UINT128_HPP
