#include <iostream>
#include <iomanip>
#include <sstream>
#include <MTToolBox/AlgorithmRecursionAndTempering.hpp>
#include <MTToolBox/AlgorithmPartialBitPattern.hpp>
//#include <UnitTest++/UnitTest++.h>
#include <NTL/GF2X.h>
#include "test_temper_searcher.hpp"

using namespace MTToolBox;
using namespace NTL;
using namespace std;

int main(void)
{
    cout << "testing recursion search and temper ...";
    Tiny32 tiny(0x8f7011ee, 0xfc78ff1f, 0, 1234);
    AlgorithmPartialBitPattern<uint32_t, 32, 1, 23, 6, false> st32;
    AlgorithmPartialBitPattern<uint32_t, 32, 1, 9, 5, true> stlsb32;
    stringstream ss;
    MersenneTwister mt;
    AlgorithmRecursionAndTempering<uint32_t> searcher(mt);
    searcher.search(tiny, st32, stlsb32, false);
    const string str = tiny.getParamString();
    if (str.find("tmat") == string::npos) {
        cout << "NG" << endl;
        return -1;
    }
    if (searcher.getWeight() >= 128) {
        cout << "NG" << endl;
        return -1;
    }
    if (searcher.getDelta() > 10) {
        cout << "NG" << endl;
        return -1;
    }
    NTL::GF2X poly = searcher.getCharacteristicPolynomial();
    if (weight(poly) != searcher.getWeight()) {
        cout << "NG" << endl;
        return -1;
    }
    cout << "ok" << endl;
    return 0;
}

#if 0
SUITE(RECURSION_SEARCH_AND_TEMPER) {
    TEST(RECURSION_SEARCH_AND_TEMPER)
    {
        Tiny32 tiny(0x8f7011ee, 0xfc78ff1f, 0, 1234);
        AlgorithmPartialBitPattern<uint32_t, uint32_t, 32, 1, 23, 6, false> st32;
        AlgorithmPartialBitPattern<uint32_t, uint32_t, 32, 1, 9, 5, true> stlsb32;
        stringstream ss;
        MersenneTwister mt;
        AlgorithmRecursionAndTempering<uint32_t> searcher(mt);
        searcher.search(tiny, st32, stlsb32, false);
        const string str = tiny.getParamString();
#if defined(DEBUG)
        cout << str << endl;
#endif
        CHECK(str.find("tmat") != string::npos);
        CHECK(searcher.getWeight() < 128);
        CHECK(searcher.getDelta() <= 10);
        NTL::GF2X poly = searcher.getCharacteristicPolynomial();
        CHECK(weight(poly) == searcher.getWeight());
    }
}
#endif
