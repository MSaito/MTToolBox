#include <iostream>
#include <iomanip>
#include <sstream>
#include <MTToolBox/AlgorithmRecursionSearch.hpp>
#include <MTToolBox/AbstractGenerator.hpp>
//#include <UnitTest++/UnitTest++.h>
#include <NTL/GF2X.h>
#include "test_generator.hpp"

using namespace MTToolBox;
using namespace NTL;
using namespace std;

int main()
{
    cout << "tesing search ...";
    Tiny32 tiny(1234);
    MersenneTwister mt;
    AlgorithmRecursionSearch<uint32_t> search(tiny, mt);
    if (search.start(10000)) {
        string str = search.getParamString();
        if (str.find("mat1") == string::npos) {
            cout << "NG" << endl;
            return -1;
        }
        GF2X poly = search.getMinPoly();
        if (deg(poly) != tiny.bitSize()) {
            cout << "NG" << endl;
            return -1;
        }
        if (search.getCount() <= 0) {
            cout << "NG" << endl;
            return -1;
        }
    } else {
        return 77;
    }
    return 0;
}
#if 0
SUITE(RECURSION_SEARCH) {
    TEST(SEARCH_U32)
    {
        Tiny32 tiny(1234);
        MersenneTwister mt;
        AlgorithmRecursionSearch<uint32_t> search(tiny, mt);
        if (search.start(10000)) {
            string str = search.getParamString();
            CHECK(str.find("mat1") != string::npos);
            GF2X poly = search.getMinPoly();
            CHECK(deg(poly) == tiny.bitSize());
            CHECK(search.getCount() > 0);
        }
    }
}
#endif
