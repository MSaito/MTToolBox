#include <iostream>
#include <iomanip>
#include <sstream>
#include <MTToolBox/AlgorithmPartialBitPattern.hpp>
//#include <UnitTest++/UnitTest++.h>
#include <NTL/GF2X.h>
#include "test_temper_searcher.hpp"

using namespace MTToolBox;
using namespace NTL;
using namespace std;

int main(void)
{
    cout << "testing temper search ...";
    Tiny32 tiny(0x8f7011ee, 0xfc78ff1f, 0, 1234);
    AlgorithmPartialBitPattern<uint32_t, uint32_t, 32, 1, 23, 6, false> st32;
    AlgorithmPartialBitPattern<uint32_t, uint32_t, 32, 1, 9, 5, true> stlsb32;
    stringstream ss;
    stlsb32(tiny, false);
    st32(tiny, false);
    string str = tiny.getParamString();
    if (str.find("tmat") == string::npos) {
        cout << "NG" << endl;
        return -1;
    }
    cout << "ok" << endl;
    return 0;
}

#if 0
SUITE(TEMPER_SEARCH) {
    TEST(TEMPER_U32)
    {
        Tiny32 tiny(0x8f7011ee, 0xfc78ff1f, 0, 1234);
        AlgorithmPartialBitPattern<uint32_t, 32, 1, 23, 6, false> st32;
        AlgorithmPartialBitPattern<uint32_t, 32, 1, 9, 5, true> stlsb32;
        stringstream ss;
        stlsb32(tiny, false);
        st32(tiny, false);
        string str = tiny.getParamString();
#if defined(DEBUG)
        cout << str << endl;
#endif
        CHECK(str.find("tmat") != string::npos);
    }
}
#endif
