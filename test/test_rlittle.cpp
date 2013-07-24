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
    bool factor = false;
    RLittle32 r(1234);
    GF2X poly;

    for (int i = 0; i < 10; i++) {
        cout << hex << r.generate() << endl;
    }
    //for (int i = 0; i < 10000000; i++) {
    for (;;) {
        minpoly(poly, r);
        if (deg(poly) > 521) {
            cout << "deg:" << dec << deg(poly) << endl;
        }
        if (deg(poly) < 521) {
            r.setUp();
            continue;
        }
        if (isIrreducible(poly)) {
            cout << "irreducible" << endl;
            r.printParam();
            irreducible = true;
        } else {
            cout << "reducible" << endl;
            r.printParam();
            reducible = true;
            if (hasFactorOfDegree(poly, 521)) {
                cout << "has a factor of degree 521" << endl;
                r.printParam();
                factor = true;
            }
        }
        if (reducible && irreducible && factor) {
            break;
        }
        r.setUp();
    }
    return 0;
}
