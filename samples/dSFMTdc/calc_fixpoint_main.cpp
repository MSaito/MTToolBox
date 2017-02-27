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
#include "calc_fixpoint.h"

using namespace MTToolBox;
using namespace std;

static bool check = false;
static bool get_fixpoint(dSFMT& dsfmt);
static bool chk_fix1(w128_t& fix, const dSFMT& dsfmt, const GF2X& irreducible,
                const GF2X& quotient);
static bool chk_fix2(const dSFMT& fix, const dSFMT& con, const GF2X& idf);

static bool chk_fix1(w128_t& fix, const dSFMT& dsfmt, const GF2X& irreducible,
                     const GF2X& quotient)
{
    GF2X a;
    GF2X b;
    GF2X d;
    dSFMT dsfmt_const(dsfmt);
    /* a*irreducible + b*quotient = d */
    XGCD(d, a, b, irreducible, quotient);
    if (deg(d) != 0) {
        cout << "failure d != 1" << endl;
        return false;
    }
#if defined(DEBUG)
    cout << "irreducible" << endl;
    printBinary(stdout, irreducible);
    cout << "quotient" << endl;
    printBinary(stdout, quotient);
    cout << "a" << endl;
    printBinary(stdout, a);
    cout << "b" << endl;
    printBinary(stdout, b);
#endif
    b *= quotient;
    GF2X idf(b);
    a *= irreducible;
    dsfmt_const.setConst();
    dSFMT dsfmt_const0(dsfmt_const);
#if defined(DEBUG)
    cout << "dsfmt_const after setConst" << endl;
    dsfmt_const.d_p();
    cout << "b" << endl;
    printBinary(stdout, b);
    cout << "dsfmt_const0" << endl;
    dsfmt_const0.d_p();
#endif
    annihilate<w128_t>(&dsfmt_const, b);
#if defined(DEBUG)
    cout << "dsfmt_const after annihilate" << endl;
    dsfmt_const.d_p();
    cout << "dsfmt_const0" << endl;
    dsfmt_const0.d_p();
#endif
    dSFMT const_L_save = dsfmt_const;
    dSFMT dsfmt_const2(dsfmt);
    dsfmt_const2.setConst();
    annihilate<w128_t>(&dsfmt_const2, a);
#if defined(DEBUG)
    cout << "dsfmt_const2 after annihilate" << endl;
    dsfmt_const2.d_p();
    cout << "dsfmt_const0" << endl;
    dsfmt_const0.d_p();
#endif
    dsfmt_const2.add(dsfmt_const);
#if defined(DEBUG)
    cout << "dsfmt_const2 after add" << endl;
    dsfmt_const2.d_p();
    cout << "dsfmt_const0" << endl;
    dsfmt_const0.d_p();
#endif
    if (!dsfmt_const2.equals(dsfmt_const0)) {
        cout << "modoranai" << endl;
        cout << "dsfmt_const2" << endl;
        dsfmt_const2.d_p();
        cout << "dsfmt_const0" << endl;
        dsfmt_const0.d_p();
        dsfmt_const2.add(dsfmt_const0);
        cout << "sa" << endl;
        dsfmt_const2.d_p();
        return false;
    }
    GF2X t1(1, 1);
    SetCoeff(t1, 0);
    /* a*irreducible + b*t1 = d */
    XGCD(d, a, b, irreducible, t1);
    if (deg(d) != 0) {
        cout << "failure d != 1" << endl;
        cout << "deg(d) = " << dec << deg(d) << endl;
        cout << "deg(irreducible) = " << dec << deg(irreducible) << endl;
        cout << "deg(t1) = " << dec << deg(t1) << endl;
        throw new logic_error("failure d != 1");
    }
    annihilate<w128_t>(&dsfmt_const, b);

    fix = dsfmt_const.getParityValue();
    if (!chk_fix2(dsfmt_const, const_L_save, idf)) {
        cout << "chk_fix2 error" << endl;
        return false;
    }
    return true;
}

static bool chk_fix2(const dSFMT& fix, const dSFMT& con, const GF2X& idf)
{
    dSFMT tmp(fix);
    tmp.generate();
    annihilate<w128_t>(&tmp, idf);
    tmp.add(&con);
    return tmp.equals(fix);
}

static bool option(dSFMT_param& param, int argc, char * argv[])
{
    int c;
    bool error = false;
    check = false;
    for (;;) {
        c = getopt(argc, argv, "hc");
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
    return true;
}

int main(int argc, char *argv[]) {
    dSFMT_param param;
    if (!option(param, argc, argv)) {
        return -1;
    }
    //cout << param.get_string() << endl;
    dSFMT dsfmt(param);
#if defined(DEBUG) && 0
    dsfmt.setConst();
    cout << "dsfmt_const after setConst" << endl;
    dsfmt.d_p();
    dsfmt.generate();
    cout << "dsfmt_const after generate" << endl;
    dsfmt.d_p();
#endif
    int mexp = dsfmt.getMexp();
    int maxdegree = dsfmt.bitSize();
    if (check) {
        printf("mexp = %d, maxdegree = %d\n", mexp, maxdegree);
    }
    w128_t w;
    w.u64[0] = 1;
    dsfmt.seed(w);
    get_fixpoint(dsfmt);
    return 0;
}

static bool get_fixpoint(dSFMT& dsfmt)
{
    GF2X poly;
    GF2X lcmpoly;
    GF2X smallpoly;
    GF2X tmp;
    GF2X rempoly;
    GF2X filler;
    w128_t fix;

    minpoly<w128_t>(poly, dsfmt);
    GF2X irreducible = poly;
    if (!hasFactorOfDegree(irreducible, dsfmt.getMexp())) {
        cout << "error does not have factor of degree "
             << dec << dsfmt.getMexp() << endl;
        cout << "degree = " << dec << deg(irreducible) << endl;
        return false;
    }
    lcmpoly = poly;
    getLCMPoly(lcmpoly, dsfmt);
    if (check) {
        printf("deg lcm poly = %ld\n", deg(lcmpoly));
    }
    DivRem(smallpoly, rempoly, lcmpoly, irreducible);
    if (deg(rempoly) >= 0) {
        cout << "deg(rempoly) = " << deg(rempoly) << endl;
        return false;
    }
    if (check) {
        printf("deg small poly = %ld\n", deg(smallpoly));
    }
    if (check) {
        if (!chk_fix1(fix, dsfmt, irreducible, smallpoly)) {
            cout << "chk_fix1 error" << endl;
            return false;
        }
    } else {
        fix = calc_fixpoint(dsfmt, irreducible, smallpoly);
    }
    dsfmt.setFixPoint(fix);
    cout << dsfmt.getParamString() << endl;
    if (check) {
        cout << "fix1 = " << setfill('0') << setw(16) << hex
             << fix.u64[0] << endl;
        cout << "fix2 = " << setfill('0') << setw(16) << hex
             << fix.u64[1] << endl;
    }
    return true;
}
