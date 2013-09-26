#ifndef XORSHIFT128_H
#define XORSHIFT128_H
#include <stdint.h>

#if defined(__cplusplus)
extern "C" {
#endif

    struct XORSHIFT128_T {
        uint32_t x;
        uint32_t y;
        uint32_t z;
        uint32_t w;
    };

    typedef struct XORSHIFT128_T xorshift128_t;

    uint32_t xorshift128_generate(xorshift128_t * xor);
    void xorshift128_init(xorshift128_t * xor, uint32_t seed);
    void xorshift128_init_by_array(xorshift128_t * xor,
                                   const uint32_t array[],
                                   int length);

#if defined(__cplusplus)
}
#endif

#endif // XORSHIFT128_H
