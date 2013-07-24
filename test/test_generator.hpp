#include <tinymt32.h>
#include <tinymt64.h>
#include <iostream>
#include <iomanip>
#include <MTToolBox/MersenneTwister.hpp>
#include <MTToolBox/abstract_generator.hpp>
#include <MTToolBox/uint128.hpp>

namespace MTToolBox {
    using namespace std;

    class Tiny32 : public U32Generator {
    public:
        Tiny32(uint32_t seed) {
            tiny.mat1 = 0x8f7011ee;
            tiny.mat2 = 0xfc78ff1f;
            tiny.tmat = 0x3793fdff;
            tinymt32_init(&tiny, seed);
        }
        uint32_t generate() {
            return tinymt32_generate_uint32(&tiny);
        }
        void seed(uint32_t value) {
            tinymt32_init(&tiny, value);
        }
        int bitSize() {
            return tinymt32_get_mexp(&tiny);
        }
    private:
        tinymt32_t tiny;
    };

    class Tiny64 : public MTToolBox::U64Generator {
    public:
        Tiny64(uint64_t seed) {
            tiny.mat1 = 0xfa051f40;
            tiny.mat2 = 0xffd0fff4;
            tiny.tmat = UINT64_C(0x58d02ffeffbfffbc);
            tinymt64_init(&tiny, seed);
        }
        uint64_t generate() {
            return tinymt64_generate_uint64(&tiny);
        }
        void seed(uint64_t value) {
            tinymt64_init(&tiny, value);
        }
        int bitSize() {
            return tinymt64_get_mexp(&tiny);
        }
    private:
        tinymt64_t tiny;
    };

    class Tiny128 : public MTToolBox::U128Generator {
    public:
        Tiny128(MTToolBox::uint128_t seed) {
            tiny.mat1 = 0xfa051f40;
            tiny.mat2 = 0xffd0fff4;
            tiny.tmat = UINT64_C(0x58d02ffeffbfffbc);
            tinymt64_init(&tiny, seed.u64[0]);
        }
        MTToolBox::uint128_t generate() {
            tinymt64_generate_uint64(&tiny);
            return MTToolBox::uint128_t(tiny.status[0], tiny.status[1]);
        }
        void seed(MTToolBox::uint128_t& value) {
            tinymt64_init(&tiny, value.u64[0]);
        }
        int bitSize() {
            return tinymt64_get_mexp(&tiny);
        }
    private:
        tinymt64_t tiny;
    };

    class RTiny32 : public MTToolBox::U32Generator {
    public:
        RTiny32(uint32_t mat1, uint32_t mat2, int sh1,
                int sh2, int sh3, int sh4, uint32_t value) {
            state = new uint32_t[length];
            setParam(mat1, mat2, sh1, sh2, sh3, sh4);
            seed(value);
        }

        RTiny32(uint32_t value) {
            state = new uint32_t[length];
            setParam(0x7c9e604d, 0xe91ac974, 2, 16, 24, 23);
            seed(value);
        }

        ~RTiny32() {
            delete[] state;
        }

        void setUp() {
            mat1 = MTToolBox::MT.next();
            mat2 = MTToolBox::MT.next();
            sh1 = MT.next() % 32;
            sh2 = MT.next() % 32;
            sh3 = MT.next() % 32;
            sh4 = MT.next() % 32;
        }

        void printParam() {
            cout << "mat1:" << hex << mat1 << endl;
            cout << "mat2:" << hex << mat2 << endl;
            cout << "sh1:" << dec << sh1 << endl;
            cout << "sh2:" << dec << sh2 << endl;
            cout << "sh3:" << dec << sh3 << endl;
            cout << "sh4:" << dec << sh4 << endl;
       }

        uint32_t generate() {
            uint32_t x;
            uint32_t y;

            x = state[0];
            y = state[1];
            y ^= y >> sh1;
            y ^= y << sh2;
            x ^= (x << sh3);
            x ^= (x >> sh4);
            x ^= y;
            state[0] = state[1];
            state[1] = state[2];
            state[2] = state[3];
            state[3] = x;
            if (x & 1) {
                state[1] ^= mat1;
                state[2] ^= mat2;
            }
            return state[0];
        }

        void seed(uint32_t value) {
            state[0] = value;
            state[1] = mat1;
            state[2] = mat2;
            state[3] = 0;
            for (int i = 1; i < min_loop; i++) {
                state[i & 3] ^= i + UINT32_C(1812433253)
                    * (state[(i - 1) & 3]
                       ^ (state[(i - 1) & 3] >> 30));
            }
        }

        int bitSize() {
            return 128;
        }
    private:
        enum {length = 4, min_loop = 8, size = 128};
        uint32_t * state;
        uint32_t mat1;
        uint32_t mat2;
        int sh1;
        int sh2;
        int sh3;
        int sh4;
        void setParam(uint32_t p_mat1, uint32_t p_mat2, int p_sh1,
                      int p_sh2, int p_sh3, int p_sh4) {
            mat1 = p_mat1;
            mat2 = p_mat2;
            sh1 = p_sh1;
            sh2 = p_sh2;
            sh3 = p_sh3;
            sh4 = p_sh4;
        }

    };
}
