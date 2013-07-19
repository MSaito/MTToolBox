#include <MTToolBox/period.hpp>
#include <MTToolBox/abstract_generator.hpp>
#include <UnitTest++/UnitTest++.h>
#include <NTL/GF2X.h>
#include <tinymt32.h>
#include <tinymt64.h>

using namespace MTToolBox;
using namespace NTL;
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

class Tiny64 : public U64Generator {
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

class Tiny128 : public U128Generator {
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
    void seed(uint128_t& value) {
        tinymt64_init(&tiny, value.u64[0]);
    }
    int bitSize() {
        return tinymt64_get_mexp(&tiny);
    }
private:
    tinymt64_t tiny;
};

TEST(PERIOD_U32)
{
    Tiny32 tiny(1234);
    GF2X poly;
    minpoly(poly, tiny);
    CHECK(deg(poly) == tiny.bitSize());
}

TEST(PERIOD_U64)
{
    Tiny64 tiny(1234);
    GF2X poly;
    minpoly(poly, tiny);
    CHECK(deg(poly) == tiny.bitSize());
}

TEST(PERIOD_U128)
{
    uint128_t seed(1234, 0);
    Tiny128 tiny(seed);
    GF2X poly;
    minpoly(poly, tiny);
    CHECK(deg(poly) == tiny.bitSize());
}

