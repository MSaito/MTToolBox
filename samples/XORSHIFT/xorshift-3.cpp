/*
 * XORSHIFT 128 の シフト量a, b, c を全件検索して最大周期となるものだけを出力する
 * あわせてv ビット精度均等分布次元のΔも出力することによって、最良のパラメータを
 * 選択できるようにする。
 */
#include <MTToolBox/EquidistributionCalculatable.hpp>
#include <MTToolBox/AlgorithmEquidistribution.hpp>
#include <MTToolBox/AlgorithmRecursionSearch.hpp>
#include <MTToolBox/AlgorithmPrimitivity.hpp>
#include <MTToolBox/Sequential.hpp>
#include <sstream>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace MTToolBox;
using namespace NTL;

class XorShift : public EquidistributionCalculatable<uint32_t> {
public:
    XorShift(uint32_t v) {
        a = 0;
        b = 0;
        c = 0;
        seed(v);
    }
    XorShift(const XorShift& src) : EquidistributionCalculatable<uint32_t>() {
        a = src.a;
        b = src.b;
        c = src.c;
        x = src.x;
        y = src.y;
        z = src.z;
        w = src.w;
    }
    uint32_t generate() {
        uint32_t t = (x ^ (x << a));
        x=y;    y=z;    z=w;
        return w = (w ^ (w >> c)) ^ (t ^ (t >> b));
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
    /*
     * この初期化はよい初期化ではない。
     * 最小多項式の計算や、均等分布次元の計算においては、
     * 全状態空間が0にならないように注意するだけでよい。
     * パラメータ決定後、実際に使用するプログラムにおいては、
     * ちゃんとした初期化を使用すること。
     */
    void seed(uint32_t v) {
        w = ~v;
        x = y = z = v;
    }
    int bitSize() const { return 128; }
    void setUpParam(ParameterGenerator& generator) {
        uint32_t r = generator.getUint32();
        a = r & 0x1f;
        r = r >> 5;
        b = r & 0x1f;
        r = r >> 5;
        c = r & 0x1f;
    }
    const std::string getHeaderString() {
        return "a, b, c";
    }
    const std::string getParamString() {
        stringstream ss;
        ss << dec << a << ",";
        ss << dec << b << ",";
        ss << dec << c;
        return ss.str();
    }
private:
    uint16_t a;
    uint16_t b;
    uint16_t c;
    uint32_t x;
    uint32_t y;
    uint32_t z;
    uint32_t w;
};

void search(Sequential<uint32_t>& seq, AlgorithmPrimitivity& ap, bool first) {
    XorShift xs(1);
    AlgorithmRecursionSearch<uint32_t> rs(xs, seq, ap);
    if (first) {
        cout << "delta:" << xs.getHeaderString() << endl;
    }
    if (rs.start(0x7fff)) {
        //cout << xs.getParamString();
        AlgorithmEquidistribution<uint32_t> eq(xs, 32);
        int veq[32];
        int delta = eq.get_all_equidist(veq);
        cout << dec << delta << ":";
        cout << xs.getParamString() << endl;
    }
}

int main() {
    Sequential<uint32_t> seq(0, 0x7fff);
    AlgorithmPrimitivity ap(prime_factors2_128_1);
    bool first = true;
    try {
        for(;;) {
            search(seq, ap, first);
            first = false;
        }
    } catch (underflow_error e) {
        return 0;
    }
    return 0;
}
