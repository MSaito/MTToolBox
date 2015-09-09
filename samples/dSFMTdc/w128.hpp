#ifndef W128_HPP
#define W128_HPP
/**
 * @file SFMTsearch.hpp
 *
 * @brief SIMD oriented Fast Mersenne Twister
 * this class is used by SFMTdc.
 *
 * This file is important. Users should not change this file,
 * except they are experts in random number generation.
 * This file is used for parameter searching.
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

#include <iostream>
#include <iomanip>
#include <MTToolBox/util.hpp>
#include <stdint.h>
#include <inttypes.h>

/**
 * @namespace sfmt
 * name space for sfmt
 */
namespace MTToolBox {
    using namespace NTL;
    using namespace std;

    union w128_t {
        uint32_t u[4];
        uint64_t u64[2];
    };
    static inline w128_t make_msb_mask(int n) {
        w128_t w;
        if (n == 128) {
            w.u64[0] = UINT64_C(0xffffffffffffffff);
            w.u64[1] = UINT64_C(0xffffffffffffffff);
        } else if (n < 64) {
            w.u64[1] = ~((~UINT64_C(0)) >> n);
            w.u64[0] = 0;
        } else {
            n = n - 64;
            w.u64[1] = UINT64_C(0xffffffffffffffff);
            w.u64[0] = ~((~UINT64_C(0)) >> n);
        }
        return w;
    }

    static inline w128_t and_mask(w128_t a, w128_t b) {
        w128_t w;
        w.u64[0] = a.u64[0] & b.u64[0];
        w.u64[1] = a.u64[1] & b.u64[1];
        return w;
    }

    template<>
    inline w128_t getOne() {
        w128_t one;
        one.u64[0] = 1;
        one.u64[1] = 0;
        return one;
    }

    template<>
    inline unsigned int getBitOfPos(w128_t bits, int pos) {
        if (pos == 0) {
            return bits.u64[0] & 1;
        } else if (pos < 64) {
            return (bits.u64[0] >> pos) & 1;
        } else if (pos == 64) {
            return bits.u64[1] & 1;
        } else {
            return (bits.u64[1] >> (pos - 64)) & 1;
        }
    }

    template<>
    inline void setBitOfPos(w128_t * bits, int pos, unsigned int b) {
        b = b & 1;
        if (pos == 0) {
            uint64_t mask = ~UINT64_C(1);
            bits->u64[0] &= mask;
            bits->u64[0] |= b;
        } else if (pos < 64) {
            uint64_t mask = ~(UINT64_C(1) << pos);
            bits->u64[0] &= mask;
            bits->u64[0] |= static_cast<uint64_t>(b) << pos;
        } else if (pos == 64) {
            uint64_t mask = ~UINT64_C(1);
            bits->u64[1] &= mask;
            bits->u64[1] |= b;
        } else {
            pos = pos - 64;
            uint64_t mask = ~(UINT64_C(1) << pos);
            bits->u64[1] &= mask;
            bits->u64[1] |= static_cast<uint64_t>(b) << pos;
        }
    }

    template<>
    inline bool isZero(w128_t x) {
        return (x.u64[0] == 0) && (x.u64[1] == 0);
    }

    template<>
    inline void setZero(w128_t& x) {
        x.u64[0] = 0;
        x.u64[1] = 0;
    }

    inline const w128_t operator>>(w128_t x, int s) {
        w128_t r;
        if (s == 0) {
            r = x;
            return r;
        } else if (s == 64) {
            r.u64[0] = x.u64[1];
            r.u64[1] = 0;
        } else if (s >= 64) {
            s = s - 64;
            r.u64[0] = x.u64[1] >> s;
            r.u64[1] = 0;
        } else {
            uint64_t w = x.u64[1] << (64 - s);
            r.u64[1] = x.u64[1] >> s;
            r.u64[0] = w | (x.u64[0] >> s);
        }
        return r;
    }

    inline const w128_t operator<<(w128_t x, int s) {
        w128_t r;
        if (s == 0) {
            r = x;
            return r;
        } else if (s == 64) {
            r.u64[1] = x.u64[0];
            r.u64[0] = 0;
        } else if (s > 64) {
            s = s - 64;
            r.u64[1] = x.u64[0] << s;
            r.u64[0] = 0;
        } else {
            uint64_t w = x.u64[0] >> (64 - s);
            r.u64[0] = x.u64[0] << s;
            r.u64[1] = w | (x.u64[1] << s);
        }
        return r;
    }

    inline const w128_t operator&(w128_t x, w128_t y) {
        w128_t r = x;
        r.u64[0] &= y.u64[0];
        r.u64[1] &= y.u64[1];
        return r;
    }

    inline const w128_t operator^(w128_t x, w128_t y) {
        w128_t r = x;
        r.u64[0] ^= y.u64[0];
        r.u64[1] ^= y.u64[1];
        return r;
    }

    inline const w128_t operator~(w128_t x) {
        w128_t r;
        r.u64[0] = ~x.u64[0];
        r.u64[1] = ~x.u64[1];
        return r;
    }

    inline w128_t& operator|=(w128_t& x, w128_t y) {
        x.u64[0] |= y.u64[0];
        x.u64[1] |= y.u64[1];
        return x;
    }

    inline w128_t& operator^=(w128_t& x, w128_t y) {
        x.u64[0] ^= y.u64[0];
        x.u64[1] ^= y.u64[1];
        return x;
    }

    inline bool operator==(const w128_t& x, const w128_t y) {
        return (x.u64[0] == y.u64[0]) && (x.u64[1] == y.u64[1]);
    }

    inline ostream& operator<<(ostream& os, w128_t x) {
        os << setfill('0');
        os << setw(16);
        os << x.u64[1];
        os << ".";
        os << setfill('0');
        os << setw(16);
        os << x.u64[0];
        return os;
    }

    static inline int calc_1pos(w128_t x)
    {
        if (isZero(x)) {
            return -1;
        }
        if (x.u64[0] == 0) {
            int64_t y = (int64_t)x.u64[1];
            y = count_bit((uint64_t)(y & -y) - 1);
            return 63 - y;
        } else {
            int64_t y = (int64_t)x.u64[0];
            y = count_bit((uint64_t)(y & -y) - 1);
            return 63 - y + 64;
        }
    }

    template<>
    inline w128_t convert(uint64_t x) {
        w128_t w;
        w.u64[0] = x;
        w.u64[1] = 0;
        return w;
    }

    template<>
    inline w128_t convert(uint32_t x) {
        w128_t w;
        w.u64[0] = x;
        w.u64[1] = 0;
        return w;
    }

    static inline w128_t reverse_bit(w128_t x) {
        w128_t w;
        w.u64[0] = reverse_bit(x.u64[1]);
        w.u64[1] = reverse_bit(x.u64[0]);
        return w;
    }
}
#endif
