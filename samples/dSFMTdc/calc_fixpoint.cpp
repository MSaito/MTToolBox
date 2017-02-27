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
#include "calc_fixpoint.h"

using namespace MTToolBox;
using namespace std;

w128_t calc_fixpoint(const dSFMT& dsfmt, const GF2X& irreducible,
                     const GF2X& quotient)
{
    GF2X a, b, d;
    dSFMT dsfmt_const(dsfmt);
    /* a*irreducible + b*quotient = d */
    XGCD(d, a, b, irreducible, quotient);
    if (deg(d) != 0) {
        cout << "failure d != 1" << endl;
        throw new logic_error("failure d != 1");
    }
    b *= quotient;
    a *= irreducible;
    dsfmt_const.setConst();
    annihilate<w128_t>(&dsfmt_const, b);

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
    return dsfmt_const.getParityValue();
}
