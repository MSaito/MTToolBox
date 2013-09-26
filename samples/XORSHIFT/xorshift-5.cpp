/*
 * xorshift-4 で求めた128ビットXORSHIFT最良の均等分布次元を再確認する。
 */
#include <MTToolBox/EquidistributionCalculatable.hpp>
#include <MTToolBox/AlgorithmEquidistribution.hpp>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace MTToolBox;
using namespace NTL;

class XorShift : public  EquidistributionCalculatable<uint32_t> {
public:
    XorShift(uint32_t v) { seed(v); }
    XorShift(const XorShift& src) : EquidistributionCalculatable<uint32_t>() {
        x = src.x;
        y = src.y;
        z = src.z;
        w = src.w;
    }
    uint32_t generate() {
        uint32_t t = x ^ (x << b);
            t ^= t >> c;
        x=y;    y=z;    z=w;
        w = (w ^ (w << a)) ^ t;
        return w;
    }
    uint32_t generate(int outBitLen) {
        uint32_t mask = 0;
        mask = (~mask) << (32 - outBitLen);
        return generate() & mask;
    }
    bool isZero() const { return x == 0 && y == 0 && z == 0 && w == 0; }
    void setZero() { x = y = z = w = 0; }
    void add(EquidistributionCalculatable<uint32_t>& other) {
        XorShift* that = dynamic_cast<XorShift *>(&other);
        if (that == 0) {
            throw std::invalid_argument(
                "the adder should have the same type as the addee.");
        }
        x ^= that->x;
        y ^= that->y;
        z ^= that->z;
        w ^= that->w;
    }
    XorShift * clone() const {
        return new XorShift(*this);
    }
    void seed(uint32_t v) {
        w = ~v;
        x = y = z = v;
    }
    int bitSize() const { return 128; }
    void setUpParam(AbstractGenerator<uint32_t>& generator){}
    const std::string getHeaderString(){return ""; }
    const std::string getParamString(){return ""; }
private:
    enum {a = 20, b = 11, c = 7};
    uint32_t x, y, z, w;
};

int main() {
    XorShift xs(1);
    AlgorithmEquidistribution<uint32_t> eq(xs, 32);
    int veq[32];
    int delta = eq.get_all_equidist(veq);
    int bitSize = xs.bitSize();
    for (int i = 0; i < 32; i++) {
        cout << "k("<< dec << setw(2) << (i + 1) << "):"
             << setw(3) << veq[i] << "  d(" << setw(2) << (i + 1) << "):"
             << setw(3) << (bitSize / (i + 1)) - veq[i] << endl;
    }
    cout << "delta:" << delta << endl;
    return 0;
}
