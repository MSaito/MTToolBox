#include <iostream>
#include <iomanip>
#include <sstream>
#include <MTToolBox/TestLinearity.hpp>
//#include <UnitTest++/UnitTest++.h>
#include "test_temper_searcher.hpp"

using namespace MTToolBox;
using namespace NTL;
using namespace std;

int main(void)
{
    cout << "testing linearity ...";
    Tiny32 tiny(0x8f7011ee, 0xfc78ff1f, 0, 1234);
    TestLinearity<uint32_t> TL;
    if (!TL(tiny)) {
        cout << "NG" << endl;
        return -1;
    }
    cout << "ok" << endl;
    return 0;
}

#if 0
SUITE(TEST_LINEARTY) {
    TEST(LINEARITY_U32)
    {
        Tiny32 tiny(0x8f7011ee, 0xfc78ff1f, 0, 1234);
        TestLinearity<uint32_t> TL;
        CHECK(TL(tiny));
    }
}
#endif
