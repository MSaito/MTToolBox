#include "xorshift128.h"

#define A 20
#define B 11
#define C 7

static void period_certification(xorshift128_t * xor);
static uint32_t ini_func1(uint32_t x);
static uint32_t ini_func2(uint32_t x);

uint32_t xorshift128_generate(xorshift128_t * xor)
{
    uint32_t t = xor->x ^ (xor->x << B);
    t ^= t >> C;
    xor->x = xor->y;
    xor->y = xor->z;
    xor->z = xor->w;
    xor->w = (xor->w ^ (xor->w << A)) ^ t;
    return xor->w;
}

void xorshift128_init(xorshift128_t * xor, uint32_t seed)
{
    uint32_t ar[4];
    ar[0] = seed;
    for (int i = 0; i < 8; i++) {
        ar[i] = UINT32_C(1812433253) * (ar[i - 1]
                                        ^ (ar[i - 1] >> 30))
            + i;
    }
    xor->x = ar[0];
    xor->y = ar[1];
    xor->z = ar[2];
    xor->w = ar[3];
}

void xorshift128_init_by_array(xorshift128_t * xor,
                               const uint32_t array[],
                               int length)
{
    const int min_loop = 8;
    const int lag = 1;
    const int mid = 1;
    const int size = 4;
    int i, j;
    int count;
    uint32_t r;
    uint32_t st[4];

    st[0] = 0;
    st[1] = 0;
    st[2] = 0;
    st[3] = 0;
    if (length + 1 > min_loop) {
        count = length + 1;
    } else {
        count = min_loop;
    }
    r = ini_func1(st[0] ^ st[mid % size]
                  ^ st[(size - 1) % size]);
    st[mid % size] += r;
    r += length;
    st[(mid + lag) % size] += r;
    st[0] = r;
    count--;
    for (i = 1, j = 0; (j < count) && (j < length); j++) {
        r = ini_func1(st[i % size]
                      ^ st[(i + mid) % size]
                      ^ st[(i + size - 1) % size]);
        st[(i + mid) % size] += r;
        r += array[j] + i;
        st[(i + mid + lag) % size] += r;
        st[i % size] = r;
        i = (i + 1) % size;
    }
    for (; j < count; j++) {
        r = ini_func1(st[i % size]
                      ^ st[(i + mid) % size]
                      ^ st[(i + size - 1) % size]);
        st[(i + mid) % size] += r;
        r += i;
        st[(i + mid + lag) % size] += r;
        st[i % size] = r;
        i = (i + 1) % size;
    }
    for (j = 0; j < size; j++) {
        r = ini_func2(st[i % size]
                      + st[(i + mid) % size]
                      + st[(i + size - 1) % size]);
        st[(i + mid) % size] ^= r;
        r -= i;
        st[(i + mid + lag) % size] ^= r;
        st[i % size] = r;
        i = (i + 1) % size;
    }
    period_certification(xor);
}

static void period_certification(xorshift128_t * xor)
{
    if (xor->x == 0 && xor->y == 0 && xor->z == 0 && xor->w == 0) {
        xor->x = ('X' << 16) | 'O';
        xor->y = ('R' << 16) | 'S';
        xor->z = ('H' << 16) | 'I';
        xor->w = ('F' << 16) | 'T';
    }
}
static uint32_t ini_func1(uint32_t x) {
    return (x ^ (x >> 27)) * UINT32_C(1664525);
}
static uint32_t ini_func2(uint32_t x) {
    return (x ^ (x >> 27)) * UINT32_C(1566083941);
}

#if defined(MAIN)
#include <stdio.h>
int main() {
    xorshift128_t xor;
    xorshift128_init(&xor, 1);
    for (int i = 0; i < 40; i++) {
        printf("%8x ", xorshift128_generate(&xor));
        if (i % 4 == 3) {
            printf("\n");
        }
    }
}
#endif
