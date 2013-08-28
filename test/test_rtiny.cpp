#include <sstream>
#include "test_generator.hpp"
#include <MTToolBox/period.hpp>

using namespace MTToolBox;
using namespace NTL;
using namespace std;

int main()
{
    bool reducible = false;
    bool irreducible = false;
    bool prime = false;
    bool factor = true;
    RTiny32 r(1234);
    GF2X poly;
    MersenneTwister mt;
    const char * table[] = {"3", "5", "17", "257", "641", "65537",
                            "274177", "6700417", "672804213107211", NULL};
#if 0
    Vec<ZZ> zz_table;
    stringstream ss;
    long length = 0;
    for (int i = 0; table[i] != NULL; i++) {
        length = i + 1;
    }
    zz_table.SetLength(length);
    for (int i = 0; table[i] != NULL; i++) {
        ss << table[i];
        ss << " ";
        ZZ w;
        cout << ss.str() << endl;
        ss >> w;
        zz_table[i] = w;
    }
#endif
    for (int i = 0; i < 10000000; i++) {
        minpoly(poly, r);
        if (deg(poly) != 128) {
            r.setUpParam(mt);
            continue;
        }
        if (isIrreducible(poly)) {
            cout << "irreducible" << endl;
            cout << r.getParamString() << endl;
            irreducible = true;
            //if (isPrime(poly, 128, zz_table)) {
            if (isPrime(poly, 128, table)) {
                cout << "prime" << endl;
                cout << r.getParamString() << endl;
            }
        } else {
            cout << "reducible" << endl;
            cout << r.getParamString() << endl;
            reducible = true;
            if (hasFactorOfDegree(poly, 127)) {
                cout << "has a factor of degree 127" << endl;
                cout << r.getParamString() << endl;
                factor = true;
            }
        }
        if (reducible && irreducible && prime && factor) {
            break;
        }
        r.setUpParam(mt);
    }
    return 0;
}
