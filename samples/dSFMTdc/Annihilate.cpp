#include "Annihilate.h"
#include "dSFMTsearch.hpp"
#include <NTL/GF2X.h>
#include <MTToolBox/period.hpp>
#include <MTToolBox/AlgorithmReducibleRecursionSearch.hpp>

using namespace NTL;
using namespace std;
using namespace MTToolBox;
static void minPolyLung(GF2X& poly, dSFMT& sf, int pos);
static void lungLCM(GF2X& poly, dSFMT& sf);

bool anni(dSFMT& sf)
{
    dSFMT gen(sf);
    GF2X poly;
    minpoly<w128_t>(poly, sf);

    GF2X irreducible = poly;
    if (!hasFactorOfDegree(irreducible, gen.getMexp())) {
        cout << "error does not have factor of degree " << dec << gen.getMexp()
             << endl;
        return false;
    }
    minpoly(poly, sf);
    if (deg(poly) != sf.bitSize()) {
        getLCMPoly(poly, sf);
    }
    //printBinary(stdout, poly);
    GF2X quotient = poly / irreducible;
#if defined(DEBUG)
    cout << "deg irreducible = " << dec << deg(irreducible) << endl;
    cout << "deg characteristic = " << dec << deg(poly)
         << endl;
    cout << "deg quotient = " << dec << deg(quotient) << endl;
#endif
    annihilate<w128_t>(&sf, quotient);
    minpoly<w128_t>(poly, sf);
    //cout << "after annihilate deg poly = " << dec << deg(poly) << endl;
    return true;
}

void getLCMPoly(GF2X& lcm, const dSFMT& sf)
{
    dSFMT gen(sf);
    lungLCM(lcm, gen);
    int bitSize = gen.bitSize();
    GF2X poly;
    for (int i = 0; i < bitSize; i++) {
        gen.setOneBit(i);
        minPolyLung(poly, gen, 0);
        LCM(lcm, lcm, poly);
        if (deg(lcm) == bitSize) {
            return;
        }
    }
}

static void lungLCM(GF2X& lcm, dSFMT& sf)
{
    GF2X poly;
    for (int i = 0; i < 128; i++) {
        minPolyLung(poly, sf, i);
        LCM(lcm, lcm, poly);
    }
}

static void minPolyLung(GF2X& poly, dSFMT& sf, int pos)
{
    Vec<GF2> v;
    int size = sf.bitSize();
    v.SetLength(2 * size);
    for (int i = 0; i < 2 * size; i++) {
        sf.generate();
        w128_t w = sf.getParityValue();
        v[i] = getBitOfPos(w, pos);
    }
    MinPolySeq(poly, v, size);
}
