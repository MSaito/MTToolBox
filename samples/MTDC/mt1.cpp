/*
 * MT19937 の v ビット精度均等分布次元を計算する。
 */
#include <MTToolBox/EquidistributionCalculatable.hpp>
#include <MTToolBox/AlgorithmEquidistribution.hpp>
#include <MTToolBox/util.hpp>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace MTToolBox;
using namespace NTL;

class MT19937 : public  EquidistributionCalculatable<uint32_t> {
public:
    MT19937(uint32_t v) { seed(v); }
    MT19937(const MT19937& src) : EquidistributionCalculatable<uint32_t>() {
        index = src.index;
        for (int i = 0; i < N; i++) {
            state[i] = src.state[i];
        }
    }
    uint32_t generate() {
        static const uint32_t matrix_a[2] = {0, UINT32_C(0x9908b0df)};
        uint32_t y;
        index = (index + 1) % N;
        y = (state[index] & upper_mask) |(state[(index + 1) % N] & lower_mask);
        y = state[(index + M) % N] ^ (y >> 1) ^ matrix_a[y & 1];
        state[index] = y;
        y ^= (y >> 11);
        y ^= (y << 7) & UINT32_C(0x9d2c5680);
        y ^= (y << 15) & UINT32_C(0xefc60000);
        y ^= (y >> 18);
        return y;
    }
    uint32_t generate(int outBitLen) {
        uint32_t mask = 0;
        mask = (~mask) << (32 - outBitLen);
        return generate() & mask;
    }

    bool isZero() const {
        if ((upper_mask & state[index]) != 0) {
            return false;
        }
        for (int i = 1; i < N; i++) {
            if (state[(index + i) % N] != 0) {
                return false;
            }
        }
        return true;
    }
    void setZero() {
        for (int i = 0; i < N; i++) {
            state[i] = 0;
        }
        index = 0;
    }
    void add(EquidistributionCalculatable<uint32_t>& other) {
        MT19937* that = dynamic_cast<MT19937 *>(&other);
        if (that == 0) {
            throw std::invalid_argument(
                "the adder should have the same type as the addee.");
        }
        for (int i = 0; i < N; i++) {
            state[(index + i) % N] ^= that->state[(that->index + i) % N];
        }
    }
    MT19937 * clone() const {
        return new MT19937(*this);
    }
    void seed(uint32_t v) {
        state[0]= v;
        for (int i = 1; i < N; i++) {
            state[i] = UINT32_C(1812433253)
                * (state[i - 1] ^ (state[i - 1] >> 30)) + i;
        }
        index = N - 1;
    }
    int bitSize() const { return 19937; }
    void setUpParam(ParameterGenerator& generator) {
        UNUSED_VARIABLE(&generator);
    }
    const std::string getHeaderString(){return ""; }
    const std::string getParamString(){return ""; }
private:
    static const uint32_t upper_mask = UINT32_C(0x80000000);
    static const uint32_t lower_mask = UINT32_C(0x7fffffff);
    enum {N = 624, M = 397};
    uint32_t state[N];
    int index;
};

int main() {
    MT19937 mt(1);
    AlgorithmEquidistribution<uint32_t> eq(mt, 32);
    int veq[32];
    int delta = eq.get_all_equidist(veq);
    int bitSize = mt.bitSize();
    for (int i = 0; i < 32; i++) {
        cout << "k("<< dec << setw(2) << (i + 1) << "):"
             << setw(5) << veq[i] << "  d(" << setw(2) << (i + 1) << "):"
             << setw(5) << (bitSize / (i + 1)) - veq[i] << endl;
    }
    cout << "delta:" << delta << endl;
    return 0;
}
