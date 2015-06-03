#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdint.h>
#include <inttypes.h>
#include <MTToolBox/AlgorithmEquidistribution.hpp>
#include <MTToolBox/AbstractGenerator.hpp>
//#include <UnitTest++/UnitTest++.h>
#include "test_generator.hpp"

using namespace MTToolBox;
using namespace std;

int main(void)
{
    cout << "testing equidistribution ...";
    Tiny32 tiny(1234);
    AlgorithmEquidistribution<uint32_t> eq(tiny, 32);
    int veq[32];
    int delta = eq.get_all_equidist(veq);
    if (delta != 0) {
        cout << "NG" << endl;
        return - 1;
    }
    cout << "ok" << endl;
    return 0;
}

#if 0
SUITE(CALC_EQUIDISTRIBUTION) {
    TEST(EQUI_U32)
    {
        Tiny32 tiny(1234);
        AlgorithmEquidistribution<uint32_t> eq(tiny, 32);
        int veq[32];
        int delta = eq.get_all_equidist(veq);
#if 0
        cout << "delta:" << delta << endl;
        int bitSize = tiny.bitSize();
        for (int i = 0; i < 32; i++) {
            cout << dec << (i + 1) << ":"
                 << veq[i] << ":"
                 << (bitSize / (i + 1)) - veq[i] << endl;
        }
#endif
        CHECK(delta == 0);
    }
}
#endif
