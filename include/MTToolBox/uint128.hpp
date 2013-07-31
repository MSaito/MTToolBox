#ifndef MTTOOLBOX_UINT128_T
#define MTTOOLBOX_UINT128_T
#include <stdint.h>

namespace MTToolBox {
    union uint128_t {
        uint32_t u32[4];
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

        uint128_t(uint32_t a, uint32_t b, uint32_t c, uint32_t d) {
            u32[0] = a;
            u32[1] = b;
            u32[2] = c;
            u32[3] = d;
        }

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
}
#endif // MTTOOLBOX_UINT128_T
