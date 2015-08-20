#include "rmt64.hpp"
#include <MTToolBox/AlgorithmEquidistribution.hpp>
#include <MTToolBox/AlgorithmReducibleRecursionSearch.hpp>
#include <MTToolBox/period.hpp>
#include <NTL/GF2X.h>
#include <errno.h>
#include <stdlib.h>
using namespace MTToolBox;
using namespace std;

int main(int argc, char * argv[])
{
    if (argc <= 1) {
        cout << argv[0] << "\n mexp,pos,mata,mskb,mskc,parity" << endl;
        return 1;
    }
    char * para = argv[1];
    int mexp = strtoull(para, &para, 10);
    para++;
    int pos = strtoull(para, &para, 10);
    para++;
    uint64_t mata = strtoull(para, &para, 16);
    para++;
    uint64_t parity = strtoull(para, &para, 16);
    para++;
    uint64_t mskb = strtoull(para, &para, 16);
    para++;
    uint64_t mskc = strtoull(para, &para, 16);

    RMT64Search rmt(mexp, pos, mata, mskb, mskc, 1234);
    rmt.setParityValue(parity);
    //cout << rmt.getParamString() << endl;
    rmt.seed(1234);
    GF2X poly;
    minpoly<uint64_t>(poly, rmt);
    //cout << "deg poly = " << dec << deg(poly) << endl;
    GF2X irreducible = poly;
    //cout << "deg irreducible = " << dec << deg(irreducible) << endl;
    //cout << "mexp = " << dec << mexp << endl;
    if (!hasFactorOfDegree(irreducible, mexp)) {
        cout << "error does not have factor of degree " << dec << mexp
             << endl;
        cout << "degree = " << dec << deg(irreducible) << endl;
        return -1;
    }
    calcCharacteristicPolynomial(&rmt, poly);
    //cout << "deg characteristic = " << dec << deg(poly) << endl;
    GF2X quotient = poly / irreducible;
    //cout << "deg irreducible = " << dec << deg(irreducible) << endl;
    //cout << "deg characteristic = " << dec << deg(poly) << endl;
    //cout << "deg quotient = " << dec << deg(quotient) << endl;
    if (!IsZero(quotient)) {
        annihilate<uint64_t>(&rmt, quotient);
    }
    minpoly<uint64_t>(poly, rmt);
    //cout << "deg poly = " << dec << deg(poly) << endl;
    int delta = 0;
    AlgorithmEquidistribution<uint64_t> re(rmt, 64, mexp);
    int veq[64];
    delta = re.get_all_equidist(veq);
    cout << rmt.getParamString();
    cout << "," << dec << delta << endl;
#if 1
    cout << "veq" << endl;
    for (int j = 0; j < 64; j++) {
        cout << "k(" << dec << (j + 1) << ") = " << dec << veq[j];
        cout << "\td(" << dec << (j + 1) << ") = " << dec
             << (mexp / (j + 1) - veq[j]) << endl;
    }
#endif
    return 0;
}
