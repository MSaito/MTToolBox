#include <sstream>
#include <MTToolBox/TestLinearity.hpp>
#include <UnitTest++.h>
#include "test_temper_searcher.hpp"

using namespace MTToolBox;
using namespace NTL;
using namespace std;

SUITE(TEST_LINEARTY) {
    TEST(LINEARITY_U32)
    {
        Tiny32 tiny(0x8f7011ee, 0xfc78ff1f, 0, 1234);
        TestLinearity<uint32_t> TL;
        CHECK(TL(tiny));
    }
}
