#include <sstream>
#include <MTToolBox/calc_equidist.hpp>
#include <MTToolBox/abstract_generator.hpp>
#include <UnitTest++/UnitTest++.h>
#include <NTL/GF2X.h>
#include "test_generator.hpp"

using namespace MTToolBox;
using namespace NTL;
using namespace std;

SUITE(CALC_EQUIDISTRIBUTION) {
    TEST(EQUI_U32)
    {
        Tiny32 tiny(1234);
        calc_equidist<Tiny32, uint32_t> eq(tiny, 32);
        int veq[32];
        int delta = eq.get_all_equidist(veq);
        int bitSize = tiny.bitSize();
#if 0
        cout << "delta:" << delta << endl;
        for (int i = 0; i < 32; i++) {
            cout << dec << (i + 1) << ":"
                 << veq[i] << ":"
                 << (bitSize / (i + 1)) - veq[i] << endl;
        }
#endif
        CHECK(delta == 0);
    }
}
