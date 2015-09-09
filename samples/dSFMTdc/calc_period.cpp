#include "dSFMTsearch.hpp"
#include "MTToolBox/period.hpp"
#include "Annihilate.h"
#include <NTL/GF2X.h>
#include <errno.h>
#include <stdlib.h>
#include <limits.h>
using namespace MTToolBox;
using namespace std;

void min_max_poly(int &min, int& max, dSFMT& sf);
int main(int argc, char * argv[])
{
    if (argc <= 1) {
        cout << argv[0] << "\n mexp,pos1,sl1,msk1,msk2,fix1,fix2,"
             << "parity1,parity2" << endl;
        return 1;
    }
    dSFMT_param params;
    char * para = argv[1];
    params.mexp = strtoul(para, &para, 10);
    para++;
    params.pos1 = strtoul(para, &para, 10);
    para++;
    params.sl1 = strtoul(para, &para, 10);
    para++;
    params.msk1 = strtoull(para, &para, 16);
    para++;
    params.msk2 = strtoull(para, &para, 16);
    para++;
    params.fix1 = strtoull(para, &para, 16);
    para++;
    params.fix2 = strtoull(para, &para, 16);
    para++;
    params.parity1 = strtoull(para, &para, 16);
    para++;
    params.parity2 = strtoull(para, &para, 16);
    dSFMT sf(params);
    cout << sf.getParamString() << endl;
    int se = 1;
    w128_t seedw;
    int min = INT_MAX;
    int max = -10;
    seedw.u[0] = se;
    sf.seed(seedw);
    min_max_poly(min, max, sf);
    cout << "berore anni ";
    cout << " min = " << dec << min;
    cout << " max = " << dec << max << endl;
    for (int i = 0; i < 2; i++) {
        for (int j = 2; j >= 1; j--) {
            GF2X poly;
            dSFMT sf2 = sf;
            sf2.setStartMode(i);
            sf2.setWeightMode(j);
            if (!anni(sf2)) {
                cout << "annihilate fail" << endl;
                return -1;
            }
            min_max_poly(min, max, sf2);
            cout << "start_mode = " << dec << i;
            cout << " weight_mode = " << dec << j;
            cout << " min = " << dec << min;
            cout << " max = " << dec << max << endl;
        }
    }
    return 0;
}

void min_max_poly(int &min, int&max, dSFMT& sf)
{
    GF2X poly;
    min = INT_MAX;
    max = -10;
    for (int i = 0; i < 52; i++) {
        minpoly<w128_t>(poly, sf, i);
        int d = deg(poly);
        if (d < min) {
            min = d;
        }
        if (d > max) {
            max = d;
        }
    }
    for (int i = 64; i < 104; i++) {
        minpoly<w128_t>(poly, sf, i);
        int d = deg(poly);
        if (d < min) {
            min = d;
        }
        if (d > max) {
            max = d;
        }
    }
}
