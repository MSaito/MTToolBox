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
    cout << rmt.getParamString() << endl;
    rmt.seed(1234);
    GF2X poly;
    minpoly<uint64_t>(poly, rmt);
    //cout << "deg poly = " << dec << deg(poly) << endl;
    GF2X irreducible = poly;
    if (!hasFactorOfDegree(irreducible, mexp)) {
        cout << "error does not have factor of degree " << dec << mexp
             << endl;
        cout << "degree = " << dec << deg(irreducible) << endl;
        return -1;
    }
    cout << "deg irreducible = " << dec << deg(irreducible) << endl;
    cout << "period >= 2^" << dec << deg(irreducible) << "-1" << endl;
    return 0;
}
