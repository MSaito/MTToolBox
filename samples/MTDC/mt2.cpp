/*
 * -*- coding:utf-8 -*-
 * MT19937 のダイナミッククリエイター
 */
#include <MTToolBox/TemperingCalculatable.hpp>
#include <MTToolBox/AlgorithmRecursionAndTempering.hpp>
#include <MTToolBox/MersenneTwister.hpp>
#include <MTToolBox/util.hpp>
#include <sstream>
#include <string>
#include <errno.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <getopt.h>
#include <MTToolBox/AlgorithmBestBits.hpp>

using namespace std;
using namespace MTToolBox;
using namespace NTL;

class MT32Search : public TemperingCalculatable<uint32_t> {
public:
    MT32Search(int mersenne_exponent, int ident, uint32_t v) {
        mexp = mersenne_exponent;
        size = mexp / 32 + 1;
        int r = size * 32 - mexp;
        upper_mask = 0;
        upper_mask = (~upper_mask) << r;
        lower_mask = ~upper_mask;
        maskb = 0;
        maskc = 0;
        id = ident;
        state = new uint32_t[size];
        reverse = false;
        seed(v);
    }
    MT32Search(int mersenne_exponent, int ident, int position,
               uint32_t matrix_a, uint32_t mask_b, uint32_t mask_c,
               uint32_t v) {
        mexp = mersenne_exponent;
        size = mexp / 32 + 1;
        int r = size * 32 - mexp;
        upper_mask = 0;
        upper_mask = (~upper_mask) << r;
        lower_mask = ~upper_mask;
        mata = matrix_a;
        pos = position;
        maskb = mask_b;
        maskc = mask_c;
        id = ident;
        state = new uint32_t[size];
        reverse = false;
        seed(v);
    }
    MT32Search(const MT32Search& src) : TemperingCalculatable<uint32_t>() {
        mexp = src.mexp;
        size = src.size;
        state = new uint32_t[size];
        index = src.index;
        reverse = src.reverse;
        id = src.id;
        pos = src.pos;
        mata = src.mata;
        maskb = src.maskb;
        maskc = src.maskc;
        upper_mask = src.upper_mask;
        lower_mask = src.lower_mask;
        for (int i = 0; i < size; i++) {
            state[i] = src.state[i];
        }
    }
    ~MT32Search() {
        delete[] state;
    }
    uint32_t generate() {
        const uint32_t matrix_a[2] = {0, mata};
        uint32_t y;
        index = (index + 1) % size;
        y = (state[index] & upper_mask)
            | (state[(index + 1) % size] & lower_mask);
        y = state[(index + pos) % size] ^ (y >> 1) ^ matrix_a[y & 1];
        state[index] = y;
        y ^= (y >> 12);
        y ^= (y << 7) & maskb;
        y ^= (y << 15) & maskc;
        y ^= (y >> 18);
        return y;
    }
    uint32_t generate(int outBitLen) {
        uint32_t u;
        if (reverse) {
            u = reverse_bit(generate());
        } else {
            u = generate();
        }
        uint32_t mask = 0;
        mask = (~mask) << (32 - outBitLen);
        return u & mask;
    }

    bool isZero() const {
        if ((upper_mask & state[index]) != 0) {
            return false;
        }
        for (int i = 1; i < size; i++) {
            if (state[(index + i) % size] != 0) {
                return false;
            }
        }
        return true;
    }

    void setZero() {
        for (int i = 0; i < size; i++) {
            state[i] = 0;
        }
        index = 0;
    }

    void add(EquidistributionCalculatable<uint32_t>& other) {
        MT32Search* that = dynamic_cast<MT32Search *>(&other);
        if (that == 0) {
            throw std::invalid_argument(
                "the adder should have the same type as the addee.");
        }
        for (int i = 0; i < size; i++) {
            state[(index + i) % size] ^= that->state[(that->index + i) % size];
        }
    }
    MT32Search * clone() const {
        return new MT32Search(*this);
    }
    void seed(uint32_t v) {
        state[0]= v;
        for (int i = 1; i < size; i++) {
            state[i] = UINT32_C(1812433253)
                * (state[i - 1] ^ (state[i - 1] >> 30)) + i;
        }
        index = size - 1;
    }

    int bitSize() const {
        return mexp;
    }

    void setUpParam(AbstractGenerator<uint32_t>& generator) {
        pos = generator.generate() % size;
        mata = generator.generate();
        if (id >= 0) {
            mata = (mata & UINT32_C(0xffff0000))
                | (id & UINT32_C(0x0000ffff));
        }
    }

    const std::string getHeaderString() {
        return "mexp, id, pos, mata, maskb, maskc";
    }

    const std::string getParamString() {
        stringstream ss;
        ss << dec << mexp << ",";
        ss << dec << id << ",";
        ss << dec << pos << ",";
        ss << hex << setw(8) << setfill('0') << mata << ",";
        ss << hex << setw(8) << setfill('0') << maskb << ",";
        ss << hex << setw(8) << setfill('0') << maskc << ",";
        return ss.str();
    }

    void setTemperingPattern(uint32_t mask, uint32_t pattern, int src_bit) {
        if (src_bit == 0) {
            maskb &= ~mask;
            maskb |= pattern & mask;
        } else {
            maskc &= ~mask;
            maskc |= pattern & mask;
        }
    }

    void setUpTempering() {
    }

    void setReverseOutput() {
        reverse = true;
    }

    void resetReverseOutput() {
        reverse = false;
    }

    bool isReverseOutput() {
        return reverse;
    }
#if defined(DEBUG)
    void debug_print() {
        cout << "mexp:" << dec << mexp << endl;
        cout << "size:" << dec << size << endl;
        cout << "reverse:" << reverse << endl;
        cout << "upper_mask:" << hex << upper_mask << endl;
        cout << "lower_mask:" << hex << upper_mask << endl;
        cout << "pos:" << dec << pos << endl;
        cout << "id:" << dec << id << endl;
        cout << "mata:" << hex << mata << endl;
        cout << "maskb:" << hex << maskb << endl;
        cout << "maskc:" << hex << maskc << endl;
        for (int i = 0; i < size; i++) {
            cout << "state[" << dec << i << "]:" << hex << state[i] << endl;
        }
        cout << "index:" << dec << index << endl;
    }
#endif
private:
    int mexp;
    int size;
    bool reverse;
    uint32_t upper_mask;
    uint32_t lower_mask;
    int pos;
    int id;
    uint32_t mata;
    uint32_t maskb;
    uint32_t maskc;
    uint32_t * state;
    int index;
};

class options {
public:
    int mexp;
    int count;
    int uid;
    int max_delta;
    int start;
    uint32_t seed;
    bool verbose;
    bool all;
    string filename;
};

void output_help(string& pgm)
{
    cout << pgm << " mexp id --count n --max max_delta --file filename" << endl;
}

bool parse_opt(options& opt, int argc, char **argv) {
    int c;
    bool error = false;
    string pgm = argv[0];
    static struct option longopts[] = {
        {"verbose", no_argument, NULL, 'v'},
        {"all", no_argument, NULL, 'a'},
        {"seed", required_argument, NULL, 'S'},
        {"file", required_argument, NULL, 'f'},
        {"start", required_argument, NULL, 's'},
        {"max", required_argument, NULL, 'm'},
        {"count", required_argument, NULL, 'c'},
        {NULL, 0, NULL, 0}};
    opt.verbose = false;
    opt.count = 1;
    opt.max_delta = 10;
    opt.all = false;
    errno = 0;
    for (;;) {
        c = getopt_long(argc, argv, "vas:m:f:c:S:", longopts, NULL);
        if (error) {
            break;
        }
        if (c == -1) {
            break;
        }
        switch (c) {
        case 'a':
            opt.all = true;
            break;
        case 'v':
            opt.verbose = true;
            break;
        case 's':
            opt.start = strtoull(optarg, NULL, 0);
            if (errno) {
                error = true;
                cerr << "start must be a number" << endl;
            }
            break;
        case 'S':
            opt.seed = strtoull(optarg, NULL, 0);
            if (errno) {
                error = true;
                cerr << "seed must be a number" << endl;
            }
            break;
        case 'f':
            opt.filename = optarg;
            break;
        case 'm':
            opt.max_delta = strtol(optarg, NULL, 10);
            if (errno) {
                error = true;
                cerr << "max must be a number" << endl;
            }
            break;
        case 'c':
            opt.count = strtoll(optarg, NULL, 10);
            if (errno) {
                error = true;
                cerr << "count must be a number" << endl;
            }
            break;
        case '?':
        default:
            error = true;
            break;
        }
    }
    argc -= optind;
    argv += optind;
    if (argc < 2) {
        error = true;
    } else {
        long mexp = strtol(argv[0], NULL, 10);
        static const int allowed_mexp[] = {521, 607, 1279, 2203,
                                           2281, 3217, 4253, 4423,
                                           9689, 9941, 11213, 19937,
                                           21701, 23209, 44497, -1};
        if (! errno) {
            bool found = false;
            for (int i = 0; allowed_mexp[i] > 0; i++) {
                if (mexp == allowed_mexp[i]) {
                    found = true;
                    break;
                }
            }
            if (! found) {
                error = true;
            }
        }
        if (errno || error){
            error = true;
            cerr << "mexp must be one of ";
            for (int i = 0; allowed_mexp[i] > 0; i++) {
                cerr << dec << allowed_mexp[i] << " ";
            }
            cerr << endl;
        }
        opt.mexp = mexp;
        long id = strtol(argv[1], NULL, 10);
        if (errno != 0 || (id < 0 || id > 0xffff)) {
            error = true;
            cerr << "id must be a number between 0 and " << 0xffff
                 << endl;
        }
        opt.uid = static_cast<uint32_t>(id);
    }
    if (!opt.filename.empty()) {
        ofstream ofs(opt.filename.c_str());
        if (ofs) {
            ofs.close();
        } else {
            error = true;
            cerr << "can't open file:" << opt.filename << endl;
        }
    }
    if (error) {
        output_help(pgm);
        return false;
    }
    return true;
}

int main(int argc, char * argv[]) {
    options opt;
    if (! parse_opt(opt, argc, argv)) {
        return -1;
    }
//    MT32Search mt(opt.mexp, opt.uid, 1);
//#aaa,mm,nn,rr,ww,wmask,umask,lmask,shift0,shift1,shiftB,shiftC,maskB,maskC
//83c30001, 8, 17, 23, 32, ffffffff, ff800000, 007fffff, 12, 18, 7, 15, 76b15f80, fdd50000
//    MT32Search mt(521, 1, 10, 0xfd6d0001, 1);
//    MT32Search mt(521, 1, 8, 0x83c30001, 0x76b15f80, 0xfdd50000, 1);
    MT32Search mt(521, 1, 8, 0x83c30001, 0, 0, 1);
#if defined(DEBUG)
    cout << mt.getParamString() << endl;
#endif
    MersenneTwister mto(opt.seed);
#if 0
    AlgorithmRecursionSearch<uint32_t> rs(mt, mto);
    if (rs.start(opt.mexp * 10)) {
        cout << mt.getHeaderString() << endl;
        cout << mt.getParamString() << endl;
    } else {
        cout << "not found" << endl;
        cout << mt.getParamString() << endl;
    }
#if defined(DEBUG)
    mt.debug_print();
#endif
#endif
#if 1
    static const int shifts[] = {7, 15};
    AlgorithmBestBits<uint32_t> tmp(32, shifts, 2, 15);
#if 0
    AlgorithmRecursionAndTempering<uint32_t> rt(mto);
// 521,1,10,fd6d0001,00000000,00000000,
//$ ./example4 32 521 1 1
//#aaa,mm,nn,rr,ww,wmask,umask,lmask,shift0,shift1,shiftB,shiftC,maskB,maskC
//83c30001, 8, 17, 23, 32, ffffffff, ff800000, 007fffff, 12, 18, 7, 15, 76b15f80, fdd50000
    if (rt.search(mt, tmp, tmp, opt.verbose, cout, true)) {
        return 0;
    }
#endif
    tmp(mt, true);
#endif
    cout << mt.getHeaderString() << endl;
    cout << mt.getParamString() << endl;
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
    return -1;
}
