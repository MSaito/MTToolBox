#include <iostream>
#include <iomanip>
#include <stdint.h>
#include <MTToolBox/MersenneTwister.hpp>
//#include <UnitTest++/UnitTest++.h>

using namespace MTToolBox;
using namespace std;

int array_seed(void);
int single_seed(void);

int main(void)
{
    cout << "testing Mersenne Twister" << endl;
    if (!array_seed()) {
        return -1;
    }
    if (!single_seed()) {
        return -1;
    }
    return 0;
}

int array_seed(void)
{
    cout << "testing array_seed ...";
    MersenneTwister mt;
    uint32_t init[4]={0x123, 0x234, 0x345, 0x456};
    int length = 4;
    mt.seed(init, length);
    uint32_t result[] = {1067595299, 955945823, 477289528,
                         4107218783, 4228976476};
    bool success = true;
    for (int i = 0; i < 5; i++) {
        if (mt.generate() != result[i]) {
            success = false;
            break;
        }
    }
    if (success) {
        cout << "ok" << endl;
        return 0;
    } else {
        cout << "NG" << endl;
        return -1;
    }
}

int single_seed(void)
{
    cout << "testing single_seed ...";
    MersenneTwister mt(5489);
    uint32_t result[] = {3499211612, 581869302, 3890346734,
                         3586334585, 545404204};
    bool success = true;
    for (int i = 0; i < 5; i++) {
        if (mt.generate() != result[i]) {
            success = false;
            break;
        }
    }
    if (success) {
        cout << "ok" << endl;
        return 0;
    } else {
        cout << "NG" << endl;
        return -1;
    }
}
#if 0
SUITE(MERSENNETWISTER) {
    TEST(ARRAY_SEED)
    {
        MersenneTwister mt;
        uint32_t init[4]={0x123, 0x234, 0x345, 0x456};
        int length = 4;
        mt.seed(init, length);
        uint32_t result[] = {1067595299, 955945823, 477289528,
                              4107218783, 4228976476};
        bool success = true;
        for (int i = 0; i < 5; i++) {
            if (mt.generate() != result[i]) {
                success = false;
                break;
            }
        }
        CHECK(success);
    }
    TEST(SINGLE_SEED)
    {
        MersenneTwister mt(5489);
        uint32_t result[] = {3499211612, 581869302, 3890346734,
                              3586334585, 545404204};
        bool success = true;
        for (int i = 0; i < 5; i++) {
            if (mt.generate() != result[i]) {
                success = false;
                break;
            }
        }
        CHECK(success);
    }
}
#endif
