#include "dSFMTsearch.hpp"
#include "Annihilate.h"
#include <errno.h>
#include <stdlib.h>
#include <getopt.h>
#include <errno.h>
#include <NTL/GF2X.h>
#include <NTL/GF2XFactoring.h>
#include <NTL/vector.h>
#include <MTToolBox/period.hpp>
#include <MTToolBox/AlgorithmReducibleRecursionSearch.hpp>
#include <MTToolBox/AlgorithmCalculateParity.hpp>
#include "printBinary.h"
#include "calc_parity.h"

using namespace MTToolBox;
using namespace std;

static bool check = false;
static bool verbose = true;

static void test_parity0(const dSFMT_param& param, GF2X& f, w128_t parity) {
    dSFMT dsfmt(param);
    GF2X mpoly;
    GF2X q;
    GF2X rem;
    Vec<GF2> vec;
    int i;
    int r;
    int result = 0;
    int count = 100;
    int mexp = dsfmt.getMexp();
    int point = 0;

    if (verbose) printf("start parity zero\n");
    dsfmt.setParityValue(parity);
    for (i = 0; i < count; i++) {
        if (verbose) printf("------\n");
        if (verbose) printf("==shoki (%d)\n", i);
        w128_t seed;
        seed.u64[0] = i + 1;
        seed.u64[1] = 0;
        dsfmt.seed(seed);
        //dsfmt.d_p();
        minpoly<w128_t>(mpoly, dsfmt);
        DivRem(q, rem, mpoly, f);
        if (deg(rem) != -1) {
            if (verbose) printf("minpoly = %ld\n", deg(mpoly));
            if (verbose) printf("rem != 0 deg rempoly = %ld\n", deg(rem));
            if (verbose) printf("deg q = %ld\n", deg(q));
            result = 0;
            point = i;
            break;
        }
        if (deg(mpoly) < mexp) {
            if (verbose) printf("mpoly = %ld\n", deg(mpoly));
        }
        if (deg(mpoly) < mexp) {
            result = 0;
            point = i;
            break;
        }
        r = dsfmt.periodCertification(true);
        if (r == 1) {
            if (verbose) printf("period certification OK\n");
        } else {
            if (verbose) printf("period certification NG -> OK 1\n");
            if (dsfmt.periodCertification(true) != 1) {
                result = 0;
                printf("period critification didn't change status 1!!\n");
                point = i;
                break;
            }
        }
        seed.u64[0] = i + 3;
        seed.u64[1] = 0;
        dsfmt.seed(seed);
        //dsfmt.d_p();
        annihilate<w128_t>(&dsfmt, f);
        if (verbose) printf("==zero\n");
        minpoly<w128_t>(mpoly, dsfmt);
        if (verbose || deg(mpoly) >= mexp) {
            printf("mpoly = %ld\n", deg(mpoly));
        }
        if (deg(mpoly) >= mexp) {
            printf("make zero state failed\n");
            result = 0;
            break;
        }
        //dsfmt.d_p();
        r = dsfmt.periodCertification(true);
        if (r == 1) {
            if (verbose) printf("period certification OK [ERROR]\n");
            dsfmt.d_p();
            result = 0;
            break;
        } else {
            if (verbose) printf("period certification NG -> OK\n");
            if (!dsfmt.periodCertification(true)) {
                result = 0;
                point = i;
                printf("period certification didn't chanege status!!\n");
                break;
            }
        }
        minpoly<w128_t>(mpoly, dsfmt);
        if (verbose || deg(mpoly) < mexp) {
            printf("mpoly = %ld\n", deg(mpoly));
        }
        if (deg(mpoly) < mexp) {
            result = 0;
            point = i;
            break;
        }
        r = dsfmt.periodCertification(true);
        if (r == 1) {
            if (verbose) printf("period certification OK\n");
        } else {
            if (verbose) printf("period certification NG -> OK\n");
            if (!dsfmt.periodCertification(true)) {
                printf("error!!\n");
                point = i;
                return;
            }
        }
        result++;
    }
    if (result) {
        printf("test successed %d / %d\n", result, count);
    } else {
        printf("test failed at count %d\n", point);
    }
}

static bool option(dSFMT_param& param, int argc, char * argv[])
{
    int c;
    bool error = false;
    check = false;
    verbose = false;
    for (;;) {
        c = getopt(argc, argv, "hcv");
        if (error) {
            break;
        }
        if (c == -1) {
            break;
        }
        switch (c) {
        case 'c':
            check = true;
            break;
        case 'v':
            verbose = true;
            break;
        case 'h':
        default:
            error = true;
            break;
        }
    }
    char * para = argv[optind];
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
    return true;
}

int main(int argc, char *argv[]) {
    dSFMT_param param;
    if (!option(param, argc, argv)) {
        return -1;
    }
    //cout << param.get_string() << endl;
    dSFMT dsfmt(param);
    int mexp = dsfmt.getMexp();
    int maxdegree = dsfmt.bitSize();
    if (check) {
        printf("mexp = %d, maxdegree = %d\n", mexp, maxdegree);
    }
    w128_t w;
    w.u64[0] = 1;
    dsfmt.seed(w);
    GF2X poly;
    GF2X smallpoly;
    GF2X tmp;
    GF2X rempoly;
    GF2X filler;
    w128_t parity;

    minpoly<w128_t>(poly, dsfmt);
    GF2X irreducible = poly;
    if (!hasFactorOfDegree(irreducible, dsfmt.getMexp())) {
        cout << "error does not have factor of degree "
             << dec << dsfmt.getMexp() << endl;
        cout << "degree = " << dec << deg(irreducible) << endl;
        return false;
    }
    parity = calc_parity(dsfmt, irreducible);
    cout << dsfmt.getParamString() << endl;
    if (check) {
        test_parity0(param, irreducible, parity);
    }
    return 0;
}
