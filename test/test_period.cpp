#include <iostream>
#include <iomanip>
#include <MTToolBox/period.hpp>
#include <MTToolBox/AbstractGenerator.hpp>
//#include <UnitTest++/UnitTest++.h>
#include <NTL/GF2X.h>
#include "test_generator.hpp"

using namespace MTToolBox;
using namespace NTL;
using namespace std;

bool period_u32(void);
bool period_u64(void);
bool is_irreducible(void);
bool is_prime(void);
bool is_prime2(void);
bool has_factor(void);

int main()
{
    cout << "testing period" << endl;
    if (!period_u32()) {
        return -1;
    }
    if (!period_u64()) {
        return -1;
    }
    if (!is_irreducible()) {
        return -1;
    }
    if (!is_prime()) {
        return -1;
    }
    if (!is_prime2()) {
        return -1;
    }
    if (!has_factor()) {
        return -1;
    }
    return 0;
}

bool period_u32(void)
{
    cout << "testing period u32 ...";
    Tiny32 tiny(1234);
    GF2X poly;
    minpoly(poly, tiny);
    if (deg(poly) == tiny.bitSize()) {
        cout << "ok" << endl;
        return true;
    }
    cout << "NG" << endl;
    return false;
}

bool period_u64(void)
{
    cout << "testing period u64 ...";
    Tiny64 tiny(1234);
    GF2X poly;
    minpoly(poly, tiny);
    if (deg(poly) == tiny.bitSize()) {
        cout << "ok" << endl;
        return true;
    }
    cout << "NG" << endl;
    return false;
}

bool is_irreducible(void)
{
    cout << "testing irreducible ...";
    Tiny64 tiny(1234);
    GF2X poly;
    minpoly(poly, tiny);
    if (!isIrreducible(poly)) {
        cout << "NG" << endl;
        return false;
    }

    // reducible
    RTiny32 rt(0x59c94057, 0xfd77d893, 15, 16, 3, 9, 1234);
    minpoly(poly, rt);
    if (isIrreducible(poly)) {
        cout << "NG" << endl;
        return false;
    }
    return true;
}

bool is_prime(void)
{
    cout << "testing prime 2 ...";
    Tiny64 tiny(1234);
    GF2X poly;
    minpoly(poly, tiny);
    if (!isPrime(poly)) {
        cout << "NG" << endl;
        return false;
    }
    cout << "ok" << endl;
    return true;
}

bool is_prime2(void)
{
    cout << "testing prime 2 ...";
    const char * table[] = {"3", "5", "17", "257", "641", "65537", "274177",
                            "6700417", "672804213107211", NULL};
    RTiny32 tiny(1234);
    GF2X poly;
    minpoly(poly, tiny);
    if (!isPrime(poly, 128, table)) {
        cout << "NG" << endl;
        return false;
    }

    // irreducible but not prime
    RTiny32 rt(0x474ba8c4, 0x3039cd1a, 31, 16, 7, 4, 1234);
    minpoly(poly, rt);
    if (isPrime(poly, 128, table)) {
        cout << "NG" << endl;
        return false;
    }
    return true;
}

bool has_factor(void)
{
    cout << "has factor ...";
    RLittle32 little(0x80903834, 7, 1, 31, 26, 26, 1234);
    GF2X poly;
    minpoly(poly, little);
    if (!hasFactorOfDegree(poly, 521)) {
        cout << "NG" << endl;
        return false;
    }

    // does not have
    RLittle32 little2(0xed6fdaa7, 7, 5, 27, 9, 12, 1234);
    minpoly(poly, little2);
    if (hasFactorOfDegree(poly, 521)) {
        cout << "NG" << endl;
        return false;
    }
    return true;
}

#if 0
SUITE(CALC_EQUIDISTRIBUTION) {
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
#if 0
    TEST(PERIOD_U128)
    {
        uint128_t seed(1234, 0);
        Tiny128 tiny(seed);
        GF2X poly;
        minpoly(poly, tiny);
        CHECK(deg(poly) == tiny.bitSize());
    }
#endif
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

    TEST(HASFACTOR)
    {
        RLittle32 little(0x80903834, 7, 1, 31, 26, 26, 1234);
        GF2X poly;
        minpoly(poly, little);
        CHECK(hasFactorOfDegree(poly, 521));

        // does not have
        RLittle32 little2(0xed6fdaa7, 7, 5, 27, 9, 12, 1234);
        minpoly(poly, little2);
        CHECK(!hasFactorOfDegree(poly, 521));
    }
}
#endif
