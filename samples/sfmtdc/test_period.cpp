#include "sfmtsearch.hpp"
#include "MTToolBox/period.hpp"
#include "Annihilate.h"
#include <NTL/GF2X.h>
#include <errno.h>
#include <stdlib.h>
#include <limits.h>
using namespace MTToolBox;
using namespace std;

void min_max_poly(int &min, int& max, sfmt& sf, int& max_weight);
int main(int argc, char * argv[])
{
    if (argc <= 1) {
        cout << argv[0] << "\n mexp,pos1,sl1,sl2,sr1,sr2,msk1,msk2,msk3,msk4,"
             << "parity1,parity2,parity3,parity4" << endl;
        return 1;
    }
    sfmt_param params;
    char * para = argv[1];
    params.mexp = static_cast<int>(strtoul(para, &para, 10));
    para++;
    params.pos1 = static_cast<int>(strtoul(para, &para, 10));
    para++;
    params.sl1 = static_cast<int>(strtoul(para, &para, 10));
    para++;
    params.sl2 = static_cast<int>(strtoul(para, &para, 10));
    para++;
    params.sr1 = static_cast<int>(strtoul(para, &para, 10));
    para++;
    params.sr2 = static_cast<int>(strtoul(para, &para, 10));
    para++;
    params.msk1 = static_cast<uint32_t>(strtoul(para, &para, 16));
    para++;
    params.msk2 = static_cast<uint32_t>(strtoul(para, &para, 16));
    para++;
    params.msk3 = static_cast<uint32_t>(strtoul(para, &para, 16));
    para++;
    params.msk4 = static_cast<uint32_t>(strtoul(para, &para, 16));
    para++;
    params.parity1 = static_cast<uint32_t>(strtoul(para, &para, 16));
    para++;
    params.parity2 = static_cast<uint32_t>(strtoul(para, &para, 16));
    para++;
    params.parity3 = static_cast<uint32_t>(strtoul(para, &para, 16));
    para++;
    params.parity4 = static_cast<uint32_t>(strtoul(para, &para, 16));
    sfmt sf(params);
    cout << sf.getParamString() << endl;
    int se = 1;
    w128_t seedw;
    int min = INT_MAX;
    int max = -10;
    int max_weight = 0;
    seedw.u[0] = static_cast<uint32_t>(se);
    sf.seed(seedw);
    min_max_poly(min, max, sf, max_weight);
    cout << "berore anni ";
    cout << " min = " << dec << min;
    cout << " max = " << dec << max << endl;
    cout << "max weight = " << dec << max_weight << endl;
    for (int i = 0; i < 4; i++) {
        for (int j = 4; j >= 1; j--) {
            GF2X poly;
            sfmt sf2 = sf;
            sf2.setStartMode(i);
            sf2.setWeightMode(j);
            if (!anni(sf2)) {
                cout << "annihilate fail" << endl;
                return -1;
            }
            min_max_poly(min, max, sf2, max_weight);
            cout << "start_mode = " << dec << i;
            cout << " weight_mode = " << dec << j;
            cout << " min = " << dec << min;
            cout << " max = " << dec << max << endl;
        }
    }
    return 0;
}

void min_max_poly(int &min, int&max, sfmt& sf, int& max_weight)
{
    GF2X poly;
    min = INT_MAX;
    max = -10;
    for (int i = 0; i < 128; i++) {
        minpoly<w128_t>(poly, sf, i);
        int d = static_cast<int>(deg(poly));
        if (d < min) {
            min = d;
        }
        if (d > max) {
            max = d;
            max_weight = static_cast<int>(weight(poly));
        }
    }
}
