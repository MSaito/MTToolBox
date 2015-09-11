#include "dSFMTsearch.hpp"
#include "Annihilate.h"
#include <errno.h>
#include <stdlib.h>
#include <getopt.h>
#include <errno.h>
#include <NTL/GF2X.h>
#include <NTL/GF2XFactoring.h>
#include <MTToolBox/period.hpp>
#include <MTToolBox/AlgorithmReducibleRecursionSearch.hpp>
#include "printBinary.h"

using namespace MTToolBox;
using namespace std;

void get_fixpoint(dSFMT& dsfmt);

static bool option(dSFMT_param& param, int argc, char * argv[]) {
    int c;
    bool error = false;
//    char *pgm = argv[0];
    for (;;) {
        c = getopt(argc, argv, "hd");
        if (error) {
            break;
        }
        if (c == -1) {
            break;
        }
        switch (c) {
        case 'd':
            //debug = true;
            break;
        case 'h':
        default:
            error = true;
            break;
        }
    }
    char * para = argv[1];
    param.mexp = strtoul(para, &para, 10);
    para++;
    param.pos1 = strtoul(para, &para, 10);
    para++;
    param.sl1 = strtoul(para, &para, 10);
    para++;
    param.msk1 = strtoull(para, &para, 16);
    para++;
    param.msk2 = strtoull(para, &para, 16);
    para++;
    param.fix1 = strtoull(para, &para, 16);
    para++;
    param.fix2 = strtoull(para, &para, 16);
    para++;
    param.parity1 = strtoull(para, &para, 16);
    para++;
    param.parity2 = strtoull(para, &para, 16);
    return true;
}

int main(int argc, char *argv[]) {
    dSFMT_param param;
    if (!option(param, argc, argv)) {
        return -1;
    }
    dSFMT dsfmt(param);
    int mexp = dsfmt.getMexp();
    int maxdegree = dsfmt.bitSize();
    cout << "succ = 1" << endl;
    cout << "mexp = " << dec << mexp;
    cout << ", maxdegree = " << dec << maxdegree << endl;
    w128_t w;
    w.u64[0] = 1;
    dsfmt.seed(w);
    cout << "pos1 = " << dec << param.pos1 << endl;
    cout << "sl1 = " << dec << param.sl1 << endl;
    cout << "msk1 = " << hex << setw(16) << setfill('0') << param.msk1 << endl;
    cout << "msk2 = " << hex << setw(16) << setfill('0') << param.msk2 << endl;
    cout << "fix1 = " << hex << setw(16) << setfill('0') << param.fix1 << endl;
    cout << "fix2 = " << hex << setw(16) << setfill('0') << param.fix2 << endl;
    cout << "pcv1 = " << hex << setw(16) << setfill('0') << param.parity1
         << endl;
    cout << "pcv2 = " << hex << setw(16) << setfill('0') << param.parity2
         << endl;
    get_fixpoint(dsfmt);
    return 0;
}

void get_fixpoint(dSFMT& dsfmt)
{
    GF2X poly;
    GF2X lcmpoly;

    minpoly<w128_t>(poly, dsfmt);
    GF2X irreducible = poly;
    if (!hasFactorOfDegree(irreducible, dsfmt.getMexp())) {
        cout << "error does not have factor of degree " << dec << dsfmt.getMexp()
             << endl;
        cout << "degree = " << dec << deg(irreducible) << endl;
        throw new logic_error("error does not have factor of degree ");
    }
    lcmpoly = poly;
    getLCMPoly(lcmpoly, dsfmt);
    printBinary(stdout, irreducible);
    printBinary(stdout, lcmpoly);
}

