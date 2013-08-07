#include <sstream>
#include <MTToolBox/recursion_search.hpp>
#include <MTToolBox/abstract_searcher.hpp>
#include <UnitTest++/UnitTest++.h>
#include <NTL/GF2X.h>
#include "test_generator.hpp"

using namespace MTToolBox;
using namespace NTL;
using namespace std;

SUITE(TEMPER_SEARCH) {
    TEST(TEMPER_U32)
    {
        Tiny32 tiny(1234);
        Search<Tiny32> s(tiny);
        stringstream ss;
        if (s.start(10000)) {
            s.printParam(ss);
            string str = ss.str();
            CHECK(str.find("mat1") != string::npos);
            GF2X poly = s.getMinPoly();
            CHECK(deg(poly) == tiny.bitSize());
            CHECK(s.getCount() > 0);
        }
    }
}
