#pragma once
#ifndef RMT64_HPP
#define RMT64_HPP

#include <MTToolBox/ReducibleTemperingCalculatable.hpp>
#include <MTToolBox/ReducibleGenerator.hpp>
#include <MTToolBox/MersenneTwister.hpp>
#include <MTToolBox/util.hpp>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace MTToolBox;
using namespace NTL;

class RMT64Search
    : virtual public ReducibleTemperingCalculatable<uint64_t> {
public:
    RMT64Search(int mersenne_exponent, uint64_t v) {
        mexp = mersenne_exponent;
        maskb = 0;
        maskc = 0;
        parity = 0;
        size = mexp / 64 + 1;
        state = new uint64_t[size];
        reverse = false;
        seed(v);
    }

    RMT64Search(int mersenne_exponent, int position,
               uint64_t matrix_a, uint64_t v) {
        mexp = mersenne_exponent;
        size = mexp / 64 + 1;
        mata = matrix_a;
        pos = position;
        maskb = 0;
        maskc = 0;
        parity = 0;
        state = new uint64_t[size];
        reverse = false;
        seed(v);
    }

    RMT64Search(int mersenne_exponent, int position,
               uint64_t matrix_a, int mb, int mc, uint64_t v) {
        mexp = mersenne_exponent;
        size = mexp / 64 + 1;
        mata = matrix_a;
        pos = position;
        maskb = mb;
        maskc = mc;
        parity = 0;
        state = new uint64_t[size];
        reverse = false;
        seed(v);
    }

    RMT64Search(const RMT64Search& src)
        : ReducibleTemperingCalculatable<uint64_t>() {
        mexp = src.mexp;
        size = src.size;
        state = new uint64_t[size];
        index = src.index;
        reverse = src.reverse;
        pos = src.pos;
        mata = src.mata;
        maskb = src.maskb;
        maskc = src.maskc;
        parity = src.parity;
        for (int i = 0; i < size; i++) {
            state[i] = src.state[i];
        }
    }

    ~RMT64Search() {
        delete[] state;
    }

    uint64_t generate() {
        const uint64_t matrix_a[2] = {0, mata};
        const uint64_t maska = UINT64_C(0x5555555555555555);
        uint64_t x;
        uint64_t y;
        index = (index + 1) % size;
        x = state[index];
        y = state[(index + pos) % size];
        y = y ^ (y << 17);
        x = y ^ (x >> 1) ^ matrix_a[(int)(x & UINT64_C(1))];
        state[index] = x;
        x ^= (x >> sh1) & maska;
        x ^= (x << sh2) & maskb;
        x ^= (x << sh3) & maskc;
        x ^= (x >> sh4);
        return x;
    }

    uint64_t generate(int outBitLen) {
        uint64_t u;
        if (reverse) {
            u = reverse_bit(generate());
        } else {
            u = generate();
        }
        uint64_t mask = 0;
        mask = (~mask) << (64 - outBitLen);
        return u & mask;
    }

    uint64_t getParityValue() const {
        return state[index % size];
    }

    void setParityValue(uint64_t parity_value) {
        state[index % size] = parity_value;
        parity = parity_value;
    }

    void setOneBit(int bitPos) {
        setZero();
        int idx = bitPos / 64;
        int r = bitPos % 64;
        if (idx >= size) {
#if defined(DEBUG)
            cout << "ERROR bitpos too big" << endl;
#endif
            return;
        }
        state[idx] = 1 << r;
#if defined(DEBUG) && 0
        debug_print();
        fflush(stdout);
#endif
    }

    bool isZero() const {
        for (int i = 0; i < size; i++) {
            if (state[i] != 0) {
                return false;
            }
        }
        return true;
    }

    void setZero() {
        for (int i = 0; i < size; i++) {
            state[i] = 0;
        }
        index = size - 1;
    }

    void add(EquidistributionCalculatable<uint64_t>& other) {
        RMT64Search* that = dynamic_cast<RMT64Search *>(&other);
        if (that == 0) {
            throw std::invalid_argument(
                "the adder should have the same type as the addee.");
        }
        for (int i = 0; i < size; i++) {
            state[(index + i) % size] ^=
                that->state[(that->index + i) % size];
        }
    }

    EquidistributionCalculatable<uint64_t> * clone() const {
        return new RMT64Search(*this);
    }

    void seed(uint64_t v) {
        state[0]= v;
        for (int i = 1; i < size; i++) {
            state[i] = i
                + UINT64_C(6364136223846793005)
                * (state[i - 1] ^ (state[i - 1] >> 62));
        }
        index = size - 1;
    }

    int bitSize() const {
        return size * 64;
    }

    int getMexp() const {
        return mexp;
    }

    void setUpParam(ParameterGenerator& generator) {
        pos = generator.getUint64() % size;
        mata = generator.getUint64();
    }

    const std::string getHeaderString() {
        return "mexp, pos, mata, parity, maskb, maskc";
    }

    const std::string getParamString() {
        maskb = maskb >> sh2;
        maskb = maskb << sh2;
        maskc = maskc >> sh3;
        maskc = maskc << sh3;
        stringstream ss;
        ss << dec << mexp << ",";
        ss << dec << pos << ",";
        ss << hex << setw(16) << setfill('0') << mata << ",";
        ss << hex << setw(16) << setfill('0') << parity << ",";
        ss << hex << setw(16) << setfill('0') << maskb << ",";
        ss << hex << setw(16) << setfill('0') << maskc;
        return ss.str();
    }

    void setTemperingPattern(uint64_t mask, uint64_t pattern, int src_bit) {
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
        cout << "pos:" << dec << pos << endl;
        cout << "mata:" << hex << mata << endl;
        cout << "parity:" << hex << parity << endl;
        cout << "maskb:" << hex << maskb << endl;
        cout << "maskc:" << hex << maskc << endl;
        for (int i = 0; i < size; i++) {
            cout << "state[" << dec << i << "]:" << hex << state[i] << endl;
        }
        cout << "index:" << dec << index << endl;
    }
#endif
private:
    enum {sh1 = 29, sh2 = 17, sh3 = 37, sh4 = 43};
    int mexp;
    int size;
    bool reverse;
    int pos;
    uint64_t mata;
    uint64_t parity;
    uint64_t maskb;
    uint64_t maskc;
    uint64_t * state;
    int index;
};

#endif // RMT_HPP
