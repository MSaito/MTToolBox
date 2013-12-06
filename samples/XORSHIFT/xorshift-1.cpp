/*
 * XORSHIFT 128
 */
#include <MTToolBox/AbstractGenerator.hpp>
#include <MTToolBox/AlgorithmPrimitivity.hpp>
#include <MTToolBox/period.hpp>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace MTToolBox;
using namespace NTL;

class XorShift : public AbstractGenerator<uint32_t> {
public:
    XorShift(uint32_t v) {
        seed(v);
    }
    uint32_t generate() {
        uint32_t t = (x ^ (x << a));
        x=y;    y=z;    z=w;
        return w = (w ^ (w >> c)) ^ (t ^ (t >> b));
    }
    void seed(uint32_t v) {
        w = ~v;
        x = y = z = v;
    }
    int bitSize() const {
        return 128;
    }
private:
    enum {a = 5, b = 14, c = 1};
    uint32_t x, y, z, w;
};

int main() {
    XorShift xs(1);
    GF2X poly;
    minpoly<uint32_t>(poly, xs);
    cout << "degree = " << deg(poly) << endl;
    if (isPrime(poly, 128, prime_factors2_128_1)) {
        cout << "period is 2^128 -1." << endl;
    } else {
        cout << "period is unknown." << endl;
    }
    return 0;
}
