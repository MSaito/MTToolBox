#include "dSFMTsearch.hpp"
#include <MTToolBox/TestLinearity.hpp>
#include <errno.h>
#include <stdlib.h>
using namespace MTToolBox;
using namespace std;
int main(int argc, char * argv[])
{
    if (argc <= 1) {
        cout << argv[0] << "\n mexp,pos1,sl1,msk1,msk2,"
             << "fix1,fix2,parity1,parity2" << endl;
        return 1;
    }
    dSFMT_param params;
    char * para = argv[1];
    params.mexp = static_cast<int>(strtoul(para, &para, 10));
    para++;
    params.pos1 = static_cast<int>(strtoul(para, &para, 10));
    para++;
    params.sl1 = static_cast<int>(strtoul(para, &para, 10));
    para++;
    params.msk1 = static_cast<uint64_t>(strtoull(para, &para, 16));
    para++;
    params.msk2 = static_cast<uint64_t>(strtoull(para, &para, 16));
    para++;
    params.fix1 = static_cast<uint64_t>(strtoull(para, &para, 16));
    para++;
    params.fix2 = static_cast<uint64_t>(strtoull(para, &para, 16));
    para++;
    params.parity1 = static_cast<uint64_t>(strtoull(para, &para, 16));
    para++;
    params.parity2 = static_cast<uint64_t>(strtoull(para, &para, 16));
    dSFMT sf(params);
    cout << sf.getParamString() << endl;
    bool success = true;
    for (int se = 1; se < 20; se++) {
        w128_t seedw;
        seedw.u[0] = static_cast<uint32_t>(se);
        sf.seed(seedw);
        for (int i = 0; i < 2; i++) {
            for (int j = 2; j >= 1; j--) {
                TestLinearity<w128_t> tl;
                dSFMT sf2 = sf;
                sf2.setStartMode(i);
                sf2.setWeightMode(j);
                if (tl(sf2)) {
                    cout << ".";
                } else {
                    cout << "Linearity Test fail" << endl;
                    success = false;
                }
            }
        }
    }
    if (success) {
        cout << "\nLinearity Test OK" << endl;
        return 0;
    }
}
