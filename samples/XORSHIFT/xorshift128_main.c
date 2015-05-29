#include "xorshift128.h"
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

