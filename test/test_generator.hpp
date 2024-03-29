#define LINEARITY_CHECK // for calc_equidist
#include <tinymt32.h>
#include <tinymt64.h>
#undef  LINEARITY_CHECK // for calc_equidist
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <stdexcept>
#include <MTToolBox/MersenneTwister.hpp>
#include <MTToolBox/AbstractGenerator.hpp>
#include <MTToolBox/EquidistributionCalculatable.hpp>
//#include <MTToolBox/uint128.hpp>

namespace MTToolBox {
    using namespace std;

    class Tiny32 : public EquidistributionCalculatable<uint32_t> {
    public:
        Tiny32(uint32_t seed) {
            tiny.mat1 = 0x8f7011ee;
            tiny.mat2 = 0xfc78ff1f;
            tiny.tmat = 0x3793fdff;
            tinymt32_init(&tiny, seed);
        }

        Tiny32(const Tiny32& that) :
            EquidistributionCalculatable<uint32_t>() {
            tiny = that.tiny;
        }

        Tiny32 * clone() const {
            return new Tiny32(*this);
        }

        uint32_t generate() {
            return tinymt32_generate_uint32(&tiny);
        }

        uint32_t generate(int outBitLen) {
            uint32_t mask = 0;
            mask = (~mask) << (32 - outBitLen);
            return tinymt32_generate_uint32(&tiny) & mask;
        }

        void seed(uint32_t value) {
            tinymt32_init(&tiny, value);
        }

        int bitSize() const {
            return tinymt32_get_mexp(&tiny);
        }

        void add(EquidistributionCalculatable<uint32_t>& other) {
            Tiny32* that = dynamic_cast<Tiny32 *>(&other);
            if(that == 0) {
                throw std::invalid_argument(
                    "the adder should have the same type as the addee.");
            }
            if (tiny.mat1 != that->tiny.mat1 ||
                tiny.mat2 != that->tiny.mat2 ||
                tiny.tmat != that->tiny.tmat) {
                throw std::invalid_argument(
                    "the adder should have the same parameter as the addee.");
            }
            for (int i = 0; i < 4; i++) {
                tiny.status[i] ^= that->tiny.status[i];
            }
        }
#if 0
        void next() {
            tinymt32_next_state(&tiny);
        }
#endif
        void setZero() {
            for (int i = 0; i < 4; i++) {
                tiny.status[i] = 0;
            }
        }

        bool isZero() const {
            if ((tiny.status[0] & TINYMT32_MASK) == 0 &&
                tiny.status[1] == 0 &&
                tiny.status[2] == 0 &&
                tiny.status[3] == 0) {
                return true;
            } else {
                return false;
            }
        }

        void setUpParam(ParameterGenerator& mt) {
            tiny.mat1 = mt.getUint32();
            tiny.mat2 = mt.getUint32();
        }

        const std::string getHeaderString() {
            return "";
        }

        const std::string getParamString() {
            std::stringstream ss;
            ss << "mat1:" << hex << tiny.mat1 << endl;
            ss << "mat2:" << hex << tiny.mat2 << endl;
            ss << "tmat:" << hex << tiny.tmat << endl;
            return ss.str();
        }

    private:
        tinymt32_t tiny;
    };

    class Tiny64 : public AbstractGenerator<uint64_t> {
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
        int bitSize() const {
            return tinymt64_get_mexp(&tiny);
        }
    private:
        tinymt64_t tiny;
    };
#if 0
    class Tiny128 : public AbstractGenerator<uint128_t> {
    public:
        Tiny128(uint128_t seed) {
            tiny.mat1 = 0xfa051f40;
            tiny.mat2 = 0xffd0fff4;
            tiny.tmat = UINT64_C(0x58d02ffeffbfffbc);
            tinymt64_init(&tiny, seed.u64[0]);
        }
        uint128_t generate() {
            tinymt64_generate_uint64(&tiny);
            return uint128_t(tiny.status[0], tiny.status[1]);
        }
        void seed(uint128_t value) {
            tinymt64_init(&tiny, value.u64[0]);
        }
        int bitSize() const {
            return tinymt64_get_mexp(&tiny);
        }
    private:
        tinymt64_t tiny;
    };
#endif
    class RTiny32 : public AbstractGenerator<uint32_t> {
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

        void setUpParam(ParameterGenerator& mt) {
            mat1 = mt.getUint32();
            mat2 = mt.getUint32();
            sh1 = mt.getUint32() % 32;
            sh2 = mt.getUint32() % 32;
            sh3 = mt.getUint32() % 32;
            sh4 = mt.getUint32() % 32;
        }

        const std::string getParamString() {
            std::stringstream ss;
            ss << "mat1:" << hex << mat1 << endl;
            ss << "mat2:" << hex << mat2 << endl;
            ss << "sh1:" << dec << sh1 << endl;
            ss << "sh2:" << dec << sh2 << endl;
            ss << "sh3:" << dec << sh3 << endl;
            ss << "sh4:" << dec << sh4 << endl;
            return ss.str();
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
                state[i & 3] ^= static_cast<uint32_t>(i)
                    + UINT32_C(1812433253)
                    * (state[(i - 1) & 3]
                       ^ (state[(i - 1) & 3] >> 30));
            }
        }

        int bitSize() const {
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

    class RLittle32 : public AbstractGenerator<uint32_t> {
    public:
        RLittle32(uint32_t mat1, int pos,
                  int sh1, int sh2, int sh3, int sh4,
                  uint32_t value) {
            state = new uint32_t[length];
            setParam(mat1, pos, sh1, sh2, sh3, sh4);
            seed(value);
        }

        RLittle32(uint32_t value) {
            state = new uint32_t[length];
            setParam(0x7c9e604d, 8, 2, 16, 24, 23);
            seed(value);
        }

        ~RLittle32() {
            delete[] state;
        }

        void setUpParam(ParameterGenerator& mt) {
            mat1 = mt.getUint32();
            pos = mt.getUint32() % length;
            sh1 = mt.getUint32() % 32;
            sh2 = mt.getUint32() % 32;
            sh3 = mt.getUint32() % 32;
            sh4 = mt.getUint32() % 32;
            seed(1234);
        }

        void printParam() {
            cout << "mat1:" << hex << mat1 << endl;
            cout << "pos:" << dec << pos << endl;
            cout << "sh1:" << dec << sh1 << endl;
            cout << "sh2:" << dec << sh2 << endl;
            cout << "sh3:" << dec << sh3 << endl;
            cout << "sh4:" << dec << sh4 << endl;
        }

        uint32_t generate() {
            uint32_t x;
            uint32_t y;

            index = (index + 1) % length;

            x = state[index % length];
            y = state[(index + pos) % length];
            y ^= y >> sh1;
            y ^= y << sh2;
            x ^= (x << sh3);
            x ^= (x >> sh4);
            x ^= y;
            state[index] = x;
            if (y & 1) {
                state[index] ^= mat1;
            }
            return state[index];
        }

        void seed(uint32_t value) {
            for (int i = 0; i < length; i++) {
                state[i] = 0;
            }
            state[0] = value;
            state[1] = mat1;
            for (int i = 1; i < length; i++) {
                state[i] ^= static_cast<uint32_t>(i)
                    + UINT32_C(1812433253)
                    * (state[i - 1]
                       ^ (state[i - 1] >> 30));
            }
            index = length - 1;
        }

        int bitSize() const {
            return size;
        }
    private:
        enum {length = 17, size = 544};
        uint32_t mat1;
        int pos;
        int sh1;
        int sh2;
        int sh3;
        int sh4;
        int index;
        uint32_t * state;
        void setParam(uint32_t p_mat1, int p_pos, int p_sh1,
                      int p_sh2, int p_sh3, int p_sh4) {
            mat1 = p_mat1;
            pos = p_pos;
            sh1 = p_sh1;
            sh2 = p_sh2;
            sh3 = p_sh3;
            sh4 = p_sh4;
        }
    };
}
