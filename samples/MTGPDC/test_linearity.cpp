#include <sstream>
#include <MTToolBox/TestLinearity.hpp>
#include <UnitTest++/UnitTest++.h>
#include "MTGP32search.hpp"
#include "mtgp_param.h"

using namespace MTToolBox;
using namespace NTL;
using namespace std;
using namespace mtgp;

SUITE(TEST_LINEARTY) {
    TEST(LINEARITY_U32)
    {
        mtgp_param<uint32_t> param;
        param.mexp = 11213;
        param.id = 0;
        param.pos = 7;
        param.sh1 = 13;
        param.sh2 = 4;
        param.mask = 0xffff8000;
        param.tbl[0] = 0x6b657038;
        param.tbl[1] = 0xbae56977;
        param.tbl[2] = 0xa3200006;
        param.tbl[3] = 0x20ac;
        mtgp32 mtgp(11213, 0);
        mtgp.set_param(param);
        TestLinearity<uint32_t> TL;
        CHECK(TL(mtgp));
    }
}
