#include <MTToolBox/period.hpp>
#include <MTToolBox/abstract_generator.hpp>
#include <UnitTest++/UnitTest++.h>
#include <NTL/GF2X.h>
#include "test_generator.hpp"

using namespace MTToolBox;
using namespace NTL;
using namespace std;

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

TEST(ISIRRECUDIBLE)
{
    Tiny64 tiny(1234);
    GF2X poly;
    minpoly(poly, tiny);
    CHECK(isIrreducible(poly));

    // reducible
    RTiny32 rt(0x59c94057, 0xfd77d893, 15, 16, 3, 9, 1234);
    minpoly(poly, rt);
    CHECK(!isIrreducible(poly));
}

TEST(ISPRIME)
{
    Tiny64 tiny(1234);
    GF2X poly;
    minpoly(poly, tiny);
    CHECK(isPrime(poly));
}

TEST(ISPRIME2)
{
    const char * table[] = {"3", "5", "17", "257", "641", "65537", "274177",
                            "6700417", "672804213107211", NULL};
    RTiny32 tiny(1234);
    GF2X poly;
    minpoly(poly, tiny);
    CHECK(isPrime(poly, 128, table));

    // irreducible but not prime
    RTiny32 rt(0x474ba8c4, 0x3039cd1a, 31, 16, 7, 4, 1234);
    minpoly(poly, rt);
    CHECK(!isPrime(poly, 128, table));
}
