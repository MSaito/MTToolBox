#include "sfmtsearch.hpp"
#include <MTToolBox/TestLinearity.hpp>
#include <errno.h>
#include <stdlib.h>
using namespace MTToolBox;
using namespace std;
int main(int argc, char * argv[])
{
    if (argc <= 1) {
        cout << argv[0] << "\n mexp,pos1,sl1,sl2,sr1,sr2,msk1,msk2,msk3,msk4,"
             << "parity1,parity2,parity3,parity4" << endl;
        return 1;
    }
    sfmt_param params;
    char * para = argv[1];
    params.mexp = strtoul(para, &para, 10);
    para++;
    params.pos1 = strtoul(para, &para, 10);
    para++;
    params.sl1 = strtoul(para, &para, 10);
    para++;
    params.sl2 = strtoul(para, &para, 10);
    para++;
    params.sr1 = strtoul(para, &para, 10);
    para++;
    params.sr2 = strtoul(para, &para, 10);
    para++;
    params.msk1 = strtoul(para, &para, 16);
    para++;
    params.msk2 = strtoul(para, &para, 16);
    para++;
    params.msk3 = strtoul(para, &para, 16);
    para++;
    params.msk4 = strtoul(para, &para, 16);
    para++;
    params.parity1 = strtoul(para, &para, 16);
    para++;
    params.parity2 = strtoul(para, &para, 16);
    para++;
    params.parity3 = strtoul(para, &para, 16);
    para++;
    params.parity4 = strtoul(para, &para, 16);
    sfmt sf(params);
    cout << sf.getParamString() << endl;
    bool success = true;
    for (int se = 1; se < 20; se++) {
        w128_t seedw;
        seedw.u[0] = se;
        sf.seed(seedw);
        for (int i = 0; i < 4; i++) {
            for (int j = 4; j >= 1; j--) {
                TestLinearity<w128_t> tl;
                sfmt sf2 = sf;
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
