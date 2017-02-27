#ifndef DSFMTSEARCH_HPP
#define DSFMTSEARCH_HPP
/**
 * @file dSFMTsearch.hpp
 *
 * @brief SIMD oriented Fast Mersenne Twister
 * this class is used by dSFMTdc.
 *
 * This file is important. Users should not change this file,
 * except they are experts in random number generation.
 * This file is used for parameter searching.
 *
 * @author Mutsuo Saito (Manieth Corp.)
 * @author Makoto Matsumoto (The University of Tokyo)
 *
 * Copyright (C) 2015 Mutsuo Saito, Makoto Matsumoto,
 * Manieth Corp. and Hiroshima University.
 * All rights reserved.
 *
 * The 3-clause BSD License is applied to this software, see
 * LICENSE.txt
 */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <sstream>
#include <MTToolBox/ReducibleGenerator.hpp>
#include <MTToolBox/MersenneTwister64.hpp>
#include <MTToolBox/util.hpp>
#include "w128.hpp"

/**
 * @namespace dSFMT
 * name space for dSFMT
 */
namespace MTToolBox {
    using namespace NTL;
    using namespace std;

    /**
     * @class dSFMT_param
     * @brief a class keeping parameters of dSFMT
     *
     * This class keeps parameters of dSFMT, and has some
     * method for outputting parameters.
     */
    class dSFMT_param {
    public:
        int mexp;
        int pos1;
        int sl1;
        uint64_t msk1;
        uint64_t msk2;
        uint64_t fix1;
        uint64_t fix2;
        uint64_t parity1;
        uint64_t parity2;

        dSFMT_param() {
            mexp = 0;
            pos1 = 0;
            sl1 = 0;
            msk1 = 0;
            msk2 = 0;
            fix1 = 0;
            fix2 = 0;
            parity1 = 0;
            parity2 = 0;
        }

        dSFMT_param(const dSFMT_param& src) {
            mexp = src.mexp;
            pos1 = src.pos1;
            sl1 = src.sl1;
            msk1 = src.msk1;
            msk2 = src.msk2;
            fix1 = src.fix1;
            fix2 = src.fix2;
            parity1 = src.parity1;
            parity2 = src.parity2;
        }
        /**
         * This method is used in output.hpp.
         * @return header line of output.
         */
        const string get_header() const {
            return "mexp, pos1, sl1, msk1, msk2, fix1, fix2, parity1, parity2";
        }

        /**
         * This method is used in output.hpp.
         * @return string of parameters
         */
        const string get_string() const {
            stringstream ss;
            ss << dec << mexp << ",";
            ss << dec << pos1 << ",";
            ss << dec << sl1 << ",";
            ss << hex << setw(16) << setfill('0') << msk1 << ",";
            ss << hex << setw(16) << setfill('0') << msk2 << ",";
            ss << hex << setw(16) << setfill('0') << fix1 << ",";
            ss << hex << setw(16) << setfill('0') << fix2 << ",";
            ss << hex << setw(16) << setfill('0') << parity1 << ",";
            ss << hex << setw(16) << setfill('0') << parity2 << ",";
            string s;
            ss >> s;
            return s;
        }

        /**
         * This method is used for DEBUG.
         * @return string of parameters.
         */
        const string get_debug_string() const {
            stringstream ss;
            ss << "mexp:" << dec << mexp << endl;
            ss << "pos1:" << dec << pos1 << endl;
            ss << "sl1:" << dec << sl1 << endl;
            ss << "msk1:" << hex << setw(16) << setfill('0') << msk1 << endl;
            ss << "msk2:" << hex << setw(16) << setfill('0') << msk2 << endl;
            ss << "fix1:" << hex << setw(16) << setfill('0') << fix1 << endl;
            ss << "fix2:" << hex << setw(16) << setfill('0') << fix2 << endl;
            ss << "parity1:" << hex << setw(16) << setfill('0') << parity1
               << endl;
            ss << "parity2:" << hex << setw(16) << setfill('0') << parity2
               << endl;
            string s;
            ss >> s;
            return s;
        }
    };

    /**
     * @class dSFMT
     * @brief DSFMT generator class used for dynamic creation
     *
     * This class is one of the main class of DSFMT dynamic creator.
     * This class is designed to be called from programs in MTToolBox,
     * but is not a subclass of some abstract class.
     * Instead, this class is passed to them as template parameters.
     */
    class dSFMT : public ReducibleGenerator<w128_t> {
    public:
        /**
         * Constructor by mexp.
         * @param mexp Mersenne Exponent
         */
        dSFMT(int mexp) {
            size = (mexp - 128) / 104 + 1;
            state = new w128_t[size];
            param.mexp = mexp;
            param.pos1 = 0;
            param.sl1 = 0;
            param.msk1 = 0;
            param.msk2 = 0;
            param.parity1 = 0;
            param.parity2 = 0;
            index = 0;
            start_mode = 0;
            weight_mode = 2;
            previous.u64[0] = 0;
            previous.u64[1] = 0;
            lung.u64[0] = 0;
            lung.u64[1] = 0;
            prefix = 0;
            fixed = false;
            fixedSL1 = 0;
        }

        ~dSFMT() {
            delete[] state;
        }

        /**
         * The copy constructor.
         * @param src The origin of copy.
         */
        dSFMT(const dSFMT& src) : param(src.param) {
            size = src.size;
            state = new w128_t[size];
            for (int i = 0; i < size; i++) {
                state[i] = src.state[i];
            }
            lung = src.lung;
            index = src.index;
            start_mode = src.start_mode;
            weight_mode = src.weight_mode;
            previous = src.previous;
            prefix = src.prefix;
            fixed = src.fixed;
            fixedSL1 = src.fixedSL1;
        }

        /**
         * Constructor by parameter.
         * @param src_param
         */
        dSFMT(const dSFMT_param& src_param) : param(src_param) {
            size = (src_param.mexp - 128) / 104 + 1;
            state = new w128_t[size];
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < 2; j++) {
                    state[i].u64[j] = 0;
                }
            }
            index = 0;
            start_mode = 0;
            weight_mode = 2;
            previous.u64[0] = 0;
            previous.u64[1] = 0;
            lung.u64[0] = 0;
            lung.u64[1] = 0;
            prefix = 0;
            fixed = false;
            fixedSL1 = 0;
        }

        EquidistributionCalculatable<w128_t> * clone() const {
            return new dSFMT(*this);
        }

        /**
         * This method initialize internal state.
         * This initialization is simple.
         * @param seed seed for initialization
         */
        void seed(w128_t seed) {
            //setZero();
            uint64_t * pstate = new uint64_t[(size + 1) * 2];
            for (int i = 0; i < (size + 1) * 2; i++) {
                pstate[i] = 0;
            }
            pstate[0] = seed.u64[0];
            pstate[1] = seed.u64[1];
            for (int i = 1; i < (size + 1) * 2; i++) {
                pstate[i] ^= i + UINT64_C(6364136223846793005)
                    * (pstate[i - 1] ^ (pstate[i - 1] >> 62));
            }
            for (int i = 0; i < size; i++) {
                state[i].u64[0] = pstate[i * 2];
                state[i].u64[1] = pstate[i * 2 + 1];
            }
            lung.u64[0] = pstate[size * 2];
            lung.u64[1] = pstate[size * 2 + 1];
            delete[] pstate;
            index = 0;
            previous.u64[0] = 0;
            previous.u64[1] = 0;
            setup_prefix();
        }

        void do_recursion(w128_t *r, w128_t *a, w128_t *b, w128_t *lung) {
            uint64_t t0;
            uint64_t t1;
            uint64_t L0;
            uint64_t L1;

            t0 = a->u64[0];
            t1 = a->u64[1];
            L0 = lung->u64[0];
            L1 = lung->u64[1];
            lung->u64[0] = (t0 << param.sl1) ^ (L1 >> 32) ^ (L1 << 32)
                ^ b->u64[0];
            lung->u64[1] = (t1 << param.sl1) ^ (L0 >> 32) ^ (L0 << 32)
                ^ b->u64[1];
            r->u64[0] = (lung->u64[0] >> sr1)
                ^ (lung->u64[0] & param.msk1) ^ t0;
            r->u64[1] = (lung->u64[1] >> sr1)
                ^ (lung->u64[1] & param.msk2) ^ t1;
        }

        /**
         * Important state transition function.
         */
        void next_state() {
            index = (index + 1) % size;
            do_recursion(&state[index],
                         &state[index],
                         &state[(index + param.pos1) % size],
                         &lung);
        }

        /**
         * Important method, generate new random number
         * @return new pseudo random number
         */
        w128_t generate() {
            next_state();
            w128_t r;
            index = index % size;
            int p = (index + size - 1) % size;
            switch (start_mode) {
            case 0:
                r.u64[0] = state[index].u64[0];
                r.u64[1] = state[index].u64[1];
                break;
            case 1:
            default:
                r.u64[0] = state[p].u64[1];
                r.u64[1] = state[index].u64[0];
                break;
            }
#if defined(DEBUG) && 0
            if (start_mode != 0) {
                cout << "start_mode = " << dec << start_mode
                     << " index = " << dec << index
                     << " p = " << dec << p << endl;
            }
#endif
            w128_t r2;
            switch (weight_mode) {
            case 2:
                r2 = r;
                break;
            case 1:
            default:
                r2.u64[0] = r.u64[0];
                r2.u64[1] = previous.u64[1];
                break;
            }
#if defined(DEBUG) && 0
            cout << "w:" << dec << weight_mode
                 << " s:" << dec << start_mode << endl;
            cout << "r:" << hex << r << endl;
            cout << "p:" << hex << previous << endl;
            cout << "2:" << hex << r2 << endl;
#endif
            previous = r;
            return r2;
        }

        void setStartMode(int mode) {
            start_mode = mode;
        }

        int getStartMode() {
            return start_mode;
        }

        void setWeightMode(int mode) {
            weight_mode = mode;
            previous.u64[0] = 0;
            previous.u64[1] = 0;
        }

        int getWeightMode() {
            return weight_mode;
        }
        /**
         * This method is called by the functions in simple_shortest_basis.hpp
         * This method returns \b bit_len bits of MSB of generated numbers
         * If reverse_bit_flag is set, bits are taken from LSB
         * @param bit_len bit length from MSB or LSB
         * @return generated numbers of bit_len
         */
        // これは呼ばれない
        w128_t generate(int bit_len) {
            w128_t w;
#if 0
            if (reverse_bit_flag) {
                w = reverse_bit(generate());
            } else {
                w = generate();
            }
#endif
            w = generate();
            w128_t mask = make_msb_mask(bit_len);
            return and_mask(w, mask);
        }

        /**
         * make parameters from given sequential number and
         * internal id
         * @param num sequential number
         */
        void setUpParam(ParameterGenerator& mt) {
            param.pos1 = mt.getUint64() % (size - 2) + 1;
            if (fixed) {
                param.sl1 = fixedSL1;
            } else {
                param.sl1 = mt.getUint64() % (52 - 1) + 1;
            }
            param.msk1 = mt.getUint64() | mt.getUint64();
            param.msk2 = mt.getUint64() | mt.getUint64();
            param.msk1 &= UINT64_C(0x000fffffffffffff);
            param.msk2 &= UINT64_C(0x000fffffffffffff);
        }

        void setZero() {
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < 2; j++) {
                    state[i].u64[j] = 0;
                }
            }
            lung.u64[0] = 0;
            lung.u64[1] = 0;
            index = 0;
            previous.u64[0] = 0;
            previous.u64[1] = 0;
        }

        /**
         * This method is called by the functions in the file
         * simple_shortest_basis.hpp
         * @return true if all elements of state is zero
         */
        bool isZero() const {
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < 2; j++) {
                    if (state[i].u64[j] != 0) {
                        return false;
                    }
                }
            }
            if ((lung.u64[0] != 0) || (lung.u64[1] != 0)) {
                return false;
            }
            switch (weight_mode) {
            case 2:
                return true;
            case 1:
            default:
                return previous.u[1] == 0;
            }
        }

        void setParityValue(w128_t parity) {
            lung = parity;
            param.parity1 = parity.u64[0];
            param.parity2 = parity.u64[1];
        }

        w128_t getParityValue() const {
            return lung;
        }

        void setOneBit(int bitPos) {
            setZero();
            if (bitPos < size * 104) {
                int idx = bitPos / 104;
                int p = (bitPos / 52) % 2;
                int r = bitPos % 52;
                state[idx].u64[p] = UINT64_C(1) << r;
            } else {
                bitPos = bitPos - size * 104;
                int p = (bitPos / 64) % 2;
                int r = bitPos % 64;
                lung.u64[p] = UINT64_C(1) << r;
            }
        }

        /**
         * This method is called by functions in the file
         * simple_shortest_basis.hpp addition of internal state as
         * GF(2) vector is possible when state transition function and
         * output function is GF(2)-linear.
         * @param that DSFMT generator added to this generator
         */
        void add(EquidistributionCalculatable<w128_t>& other) {
            dSFMT *that = dynamic_cast<dSFMT *>(&other);
            if (that == 0) {
                throw std::invalid_argument(
                    "the adder should have same type as the addee.");
            }
            this->add(that);
        }

        void add(const dSFMT * that) {
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < 2; j++) {
                    state[(i + index) % size].u64[j]
                        ^= that->state[(i + that->index) % size].u64[j];
                }
            }
            lung ^= that->lung;
            previous ^= that->previous;
        }

        int getMexp() const {
            return param.mexp;
        }

        int bitSize() const {
            return size * 104 + 128;
        }

        const std::string getHeaderString() {
            return param.get_header();
        }

        const std::string getParamString() {
            return param.get_string();
        }

        /**
         * This method is called by the functions in search_temper.hpp
         * to calculate the equidistribution properties from LSB
         */
        void set_reverse_bit() {
            //reverse_bit_flag = true;
        }

        /**
         * This method is called by the functions in search_temper.hpp
         * to reset the reverse_bit_flag
         */
        void reset_reverse_bit() {
            //reverse_bit_flag = false;
        }

        void setPrefix(uint64_t p) {
            prefix = p;
        }

        void setConst() {
            const uint64_t high = UINT64_C(0x3ff0000000000000);
            dSFMT tmp(*this);
            tmp.setZero();
            tmp.setPrefix(high);
            tmp.setup_prefix();
            tmp.next_state();
            setZero();
            setPrefix(high);
            setup_prefix();
            add(tmp);
            setPrefix(0);
        }

        bool equals(const dSFMT& that) {
            if ((lung.u64[0] != that.lung.u64[0])
                || (lung.u64[1] != that.lung.u64[1])) {
                return false;
            }
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < 2; j++) {
                    if (state[(index + i) % size].u64[j]
                        != that.state[(that.index + i) % size].u64[j]) {
                        return false;
                    }
                }
            }
            return true;
        }

        void setFixPoint(const w128_t& fix) {
            param.fix1 = fix.u64[0];
            param.fix2 = fix.u64[1];
        }

        void d_p() {
            cout << "index = " << dec << index << endl;
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < 2; j++) {
                    cout << setfill('0') << setw(16) << hex << state[i].u64[j]
                         << " ";
                }
                cout << endl;
            }
            for (int j = 0; j < 2; j++) {
                cout << setfill('0') << setw(16) << hex << lung.u64[j]
                     << " ";
            }
            cout << endl;
        }

        int periodCertification(bool noFix = false) {
            w128_t tmp;
            w128_t parity;
            parity.u64[0] = param.parity1;
            parity.u64[1] = param.parity2;
            w128_t fix;
            fix.u64[0] = param.fix1;
            fix.u64[1] = param.fix2;
            if (noFix) {
                tmp = lung & parity;
            } else {
                tmp = lung ^ fix;
                tmp = tmp & parity;
            }
            int c = count_bit(tmp.u64[0]);
            c += count_bit(tmp.u64[1]);
            c &= 1;
            if (c == 1) {
                return 1;
            }
            if ((parity.u64[1] & 1) == 1) {
                lung.u64[1] ^= 1;
                return 0;
            }
            for (int i = 1; i >= 0; i--) {
                uint64_t work = 1;
                for (int j = 0; j < 64; j++) {
                    if ((work & parity.u64[i]) != 0) {
                        lung.u64[i] ^= work;
                        return 0;
                    }
                    work = work << 1;
                }
            }
            return 0;
        }
        void setFixed(bool value) {
            fixed = value;
        }
        void setFixedSL1(int value) {
            fixedSL1 = value;
        }
    private:
        dSFMT& operator=(const dSFMT&) {
            throw std::logic_error("can't assign");
        }
        void setup_prefix() {
            const uint64_t clear = UINT64_C(0x000fffffffffffff);
            for (int i = 0; i < size; i++) {
                state[i].u64[0] &= clear;
                state[i].u64[1] &= clear;
            }
            if (prefix != 0) {
                for (int i = 0; i < size; i++) {
                    state[i].u64[0] |= prefix;
                    state[i].u64[1] |= prefix;
                }
            }
        }
        enum {sr1 = 12};
        int fixedSL1;
        bool fixed;
        int size;
        int index;
        int start_mode;
        int weight_mode;
        w128_t * state;
        w128_t lung;
        dSFMT_param param;
        w128_t previous;
        uint64_t prefix;
    };
}

#endif
