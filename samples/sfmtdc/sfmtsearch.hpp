#ifndef SFMTSEARCH_HPP
#define SFMTSEARCH_HPP
/**
 * @file SFMTsearch.hpp
 *
 * @brief SIMD oriented Fast Mersenne Twister
 * this class is used by SFMTdc.
 *
 * This file is important. Users should not change this file,
 * except they are experts in random number generation.
 * This file is used for parameter searching.
 *
 * @author Mutsuo Saito (Manieth Copr.)
 * @author Makoto Matsumoto (The University of Tokyo)
 *
 * Copyright (C) 2015 Mutsuo Saito, Makoto Matsumoto,
 * Manieth Copr. and Hiroshima University.
 * All rights reserved.
 *
 * The 3-clause BSD License is applied to this software, see
 * LICENSE.txt
 */

#include "w128.hpp"
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <sstream>
#include <MTToolBox/ReducibleGenerator.hpp>
#include <MTToolBox/MersenneTwister.hpp>
#include <MTToolBox/util.hpp>

/**
 * @namespace sfmt
 * name space for sfmt
 */
namespace MTToolBox {
//    using namespace MTToolBox;
    using namespace NTL;
    using namespace std;

    /**
     * @class sfmt_param
     * @brief a class keeping parameters of sfmt
     *
     * This class keeps parameters of sfmt, and has some
     * method for outputting parameters.
     */
    class sfmt_param {
    public:
        int mexp;
        int pos1;
        int sl1;
        int sl2;
        int sr1;
        int sr2;
        uint32_t msk1;
        uint32_t msk2;
        uint32_t msk3;
        uint32_t msk4;
        uint32_t parity1;
        uint32_t parity2;
        uint32_t parity3;
        uint32_t parity4;

        sfmt_param() {
            mexp = 0;
            pos1 = 0;
            sl1 = 0;
            sl2 = 0;
            sr1 = 0;
            sr2 = 0;
            msk1 = 0;
            msk2 = 0;
            msk3 = 0;
            msk4 = 0;
            parity1 = 0;
            parity2 = 0;
            parity3 = 0;
            parity4 = 0;
        }

        sfmt_param(const sfmt_param& src) {
            mexp = src.mexp;
            pos1 = src.pos1;
            sl1 = src.sl1;
            sl2 = src.sl2;
            sr1 = src.sr1;
            sr2 = src.sr2;
            msk1 = src.msk1;
            msk2 = src.msk2;
            msk3 = src.msk3;
            msk4 = src.msk4;
            parity1 = src.parity1;
            parity2 = src.parity2;
            parity3 = src.parity3;
            parity4 = src.parity4;
        }
        /**
         * This method is used in output.hpp.
         * @return header line of output.
         */
        const string get_header() const {
            return "mexp, pos1, sl1, sl2, sr1, sr2, msk1, msk2, msk3, msk4"
                ", parity1, parity2, parity3, parity4";
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
            ss << dec << sl2 << ",";
            ss << dec << sr1 << ",";
            ss << dec << sr2 << ",";
            ss << hex << setw(8) << setfill('0') << msk1 << ",";
            ss << hex << setw(8) << setfill('0') << msk2 << ",";
            ss << hex << setw(8) << setfill('0') << msk3 << ",";
            ss << hex << setw(8) << setfill('0') << msk4 << ",";
            ss << hex << setw(8) << setfill('0') << parity1 << ",";
            ss << hex << setw(8) << setfill('0') << parity2 << ",";
            ss << hex << setw(8) << setfill('0') << parity3 << ",";
            ss << hex << setw(8) << setfill('0') << parity4 << ",";
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
            ss << "sl2:" << dec << sl2 << endl;
            ss << "sr1:" << dec << sr1 << endl;
            ss << "sr2:" << dec << sr2 << endl;
            ss << "msk1:" << hex << setw(8) << setfill('0') << msk1 << endl;
            ss << "msk2:" << hex << setw(8) << setfill('0') << msk2 << endl;
            ss << "msk3:" << hex << setw(8) << setfill('0') << msk3 << endl;
            ss << "msk4:" << hex << setw(8) << setfill('0') << msk4 << endl;
            string s;
            ss >> s;
            return s;
        }
    };

    /**
     * @class sfmt
     * @brief SFMT generator class used for dynamic creation
     *
     * This class is one of the main class of SFMT dynamic creator.
     * This class is designed to be called from programs in MTToolBox,
     * but is not a subclass of some abstract class.
     * Instead, this class is passed to them as template parameters.
     */
    class sfmt : public ReducibleGenerator<w128_t> {
    public:
        /**
         * Constructor by mexp.
         * @param mexp Mersenne Exponent
         */
        sfmt(int mexp) {
            size = mexp / 128 + 1;
            state = new w128_t[size];
            param.mexp = mexp;
            param.pos1 = 0;
            param.sl1 = 0;
            param.sl2 = 0;
            param.sr1 = 0;
            param.sr2 = 0;
            param.msk1 = 0;
            param.msk2 = 0;
            param.msk3 = 0;
            param.msk4 = 0;
            param.parity1 = 0;
            param.parity2 = 0;
            param.parity3 = 0;
            param.parity4 = 0;
            index = 0;
            reverse_bit_flag = false;
            start_mode = 0;
            weight_mode = 4;
            previous.u64[0] = 0;
            previous.u64[1] = 0;
        }

        ~sfmt() {
            delete[] state;
        }

        /**
         * The copy constructor.
         * @param src The origin of copy.
         */
        sfmt(const sfmt& src) : param(src.param) {
            size = src.size;
            state = new w128_t[size];
            for (int i = 0; i < size; i++) {
                state[i] = src.state[i];
            }
            index = src.index;
            start_mode = src.start_mode;
            weight_mode = src.weight_mode;
            reverse_bit_flag = src.reverse_bit_flag;
            previous = src.previous;
        }

        /**
         * Constructor by parameter.
         * @param src_param
         */
        sfmt(const sfmt_param& src_param) : param(src_param) {
            size = src_param.mexp / 128 + 1;
            state = new w128_t[size];
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < 2; j++) {
                    state[i].u64[j] = 0;
                }
            }
            index = 0;
            start_mode = 0;
            weight_mode = 4;
            previous.u64[0] = 0;
            previous.u64[1] = 0;
            reverse_bit_flag = false;
        }

        EquidistributionCalculatable<w128_t> * clone() const {
            return new sfmt(*this);
        }

        /**
         * This method initialize internal state.
         * This initialization is simple.
         * @param seed seed for initialization
         */
        void seed(w128_t seed) {
            setZero();
            state[0] = seed;
            uint32_t * pstate = &state[0].u[0];
            for (int i = 1; i < size * 4; i++) {
                pstate[i] ^= i + UINT32_C(1812433253)
                    * (pstate[i - 1] ^ (pstate[i - 1] >> 30));
            }
            index = 0;
            previous.u64[0] = 0;
            previous.u64[1] = 0;
        }

        void rshift128(w128_t *out, w128_t const *in, int shift) {
            uint64_t th, tl, oh, ol;

            th = ((uint64_t)in->u[3] << 32) | ((uint64_t)in->u[2]);
            tl = ((uint64_t)in->u[1] << 32) | ((uint64_t)in->u[0]);

            oh = th >> (shift * 8);
            ol = tl >> (shift * 8);
            ol |= th << (64 - shift * 8);
            out->u[1] = (uint32_t)(ol >> 32);
            out->u[0] = (uint32_t)ol;
            out->u[3] = (uint32_t)(oh >> 32);
            out->u[2] = (uint32_t)oh;
        }

        void lshift128(w128_t *out, w128_t const *in, int shift) {
            uint64_t th, tl, oh, ol;

            th = ((uint64_t)in->u[3] << 32) | ((uint64_t)in->u[2]);
            tl = ((uint64_t)in->u[1] << 32) | ((uint64_t)in->u[0]);

            oh = th << (shift * 8);
            ol = tl << (shift * 8);
            oh |= tl >> (64 - shift * 8);
            out->u[1] = (uint32_t)(ol >> 32);
            out->u[0] = (uint32_t)ol;
            out->u[3] = (uint32_t)(oh >> 32);
            out->u[2] = (uint32_t)oh;
        }

        void do_recursion(w128_t *r, w128_t *a, w128_t *b,
                          w128_t *c, w128_t *d) {
            w128_t x;
            w128_t y;

            lshift128(&x, a, param.sl2);
            rshift128(&y, c, param.sr2);
            r->u[0] = a->u[0] ^ x.u[0] ^ ((b->u[0] >> param.sr1) & param.msk1)
                ^ y.u[0] ^ (d->u[0] << param.sl1);
            r->u[1] = a->u[1] ^ x.u[1] ^ ((b->u[1] >> param.sr1) & param.msk2)
                ^ y.u[1] ^ (d->u[1] << param.sl1);
            r->u[2] = a->u[2] ^ x.u[2] ^ ((b->u[2] >> param.sr1) & param.msk3)
                ^ y.u[2] ^ (d->u[2] << param.sl1);
            r->u[3] = a->u[3] ^ x.u[3] ^ ((b->u[3] >> param.sr1) & param.msk4)
                ^ y.u[3] ^ (d->u[3] << param.sl1);
        }

        /**
         * Important state transition function.
         */
        void next_state() {
            index = (index + 1) % size;
            do_recursion(&state[index],
                         &state[index],
                         &state[(index + param.pos1) % size],
                         &state[(index + size - 2) % size],
                         &state[(index + size - 1) % size]);
            //index = (index + 1) % size;
        }

        /**
         * Important method, generate new random number
         * @return new pseudo random number
         */
        w128_t generate() {
#if defined(DEBUG) && 0
            w128_t s[size];
            for (int i = 0; i < size; i++) {
                s[i] = state[i];
            }
#endif
            next_state();
            w128_t r;
            index = index % size;
            int p = (index + size - 1) % size;
            switch (start_mode) {
            case 0:
                r.u[0] = state[index].u[0];
                r.u[1] = state[index].u[1];
                r.u[2] = state[index].u[2];
                r.u[3] = state[index].u[3];
                break;
            case 1:
                r.u[0] = state[p].u[1];
                r.u[1] = state[p].u[2];
                r.u[2] = state[p].u[3];
                r.u[3] = state[index].u[0];
                break;
            case 2:
                r.u[0] = state[p].u[2];
                r.u[1] = state[p].u[3];
                r.u[2] = state[index].u[0];
                r.u[3] = state[index].u[1];
                break;
            case 3:
            default:
                r.u[0] = state[p].u[3];
                r.u[1] = state[index].u[0];
                r.u[2] = state[index].u[1];
                r.u[3] = state[index].u[2];
                break;
            }
#if defined(DEBUG) && 0
            if (start_mode != 0) {
                cout << "start_mode = " << dec << start_mode
                     << " index = " << dec << index
                     << " p = " << dec << p << endl;
            }
#endif
#if defined(DEBUG) && 0
            if ((r.u64[1] == 0 || r.u64[0] == 0) &&
                (r.u64[1] != 0 || r.u64[0] != 0)) {
                cout << "sfmt generate " << hex << r << endl;
                cout << "start_mode = " << dec << start_mode << endl;
                cout << "index = " << dec << index << endl;
                cout << "param " << getParamString() << endl;
                for (int i = 0; i < size; i++) {
                    cout << "s[" << dec << i << "] = " << hex << s[i]
                         << endl;
                }
                for (int i = 0; i < size; i++) {
                    cout << "state[" << dec << i << "] = " << hex << state[i]
                         << endl;
                }
            }
#endif
            w128_t r2;
            switch (weight_mode) {
            case 4:
                r2 = r;
                break;
            case 3:
                r2.u[0] = r.u[0];
                r2.u[1] = r.u[1];
                r2.u[2] = r.u[2];
                r2.u[3] = previous.u[3];
                break;
            case 2:
                r2.u[0] = r.u[0];
                r2.u[1] = r.u[1];
                r2.u[2] = previous.u[2];
                r2.u[3] = previous.u[3];
                break;
            case 1:
            default:
                r2.u[0] = r.u[0];
                r2.u[1] = previous.u[1];
                r2.u[2] = previous.u[2];
                r2.u[3] = previous.u[3];
                break;
            }
#if defined(DEBUG)
            cout << "w:" << dec << weight_mode
                 << " s:" << dec << start_mode << endl;
            cout << "r:" << hex << r << endl;
            cout << "p:" << hex << previous << endl;
            cout << "2:" << hex << r2 << endl;
#endif
            previous = r;
//            next_state();
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
        w128_t generate(int bit_len) {
            w128_t w;
            if (reverse_bit_flag) {
                w = reverse_bit(generate());
            } else {
                w = generate();
            }
            w128_t mask = make_msb_mask(bit_len);
            return and_mask(w, mask);
        }

        /**
         * make parameters from given sequential number and
         * internal id
         * @param num sequential number
         */
        void setUpParam(ParameterGenerator& mt) {
            param.pos1 = mt.getUint32() % (size - 2) + 1;
#if defined(SFMT_PARAM_FIXED)
            // These parameters are not best ones.
            param.sl1 = 19;
            param.sl2 = 1;
            param.sr1 = 7;
            param.sr2 = 1;
#else
            param.sl1 = mt.getUint32() % (32 - 1) + 1;
            param.sl2 = (mt.getUint32() % 4) * 2 + 1;
            param.sr1 = mt.getUint32() % (32 - 1) + 1;
            param.sr2 = (mt.getUint32() % 4) * 2 + 1;
#endif
            param.msk1 = mt.getUint32() | mt.getUint32();
            param.msk2 = mt.getUint32() | mt.getUint32();
            param.msk3 = mt.getUint32() | mt.getUint32();
            param.msk4 = mt.getUint32() | mt.getUint32();
        }

        void setZero() {
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < 2; j++) {
                    state[i].u64[j] = 0;
                }
            }
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
            switch (weight_mode) {
            case 4:
                return true;
            case 3:
                return previous.u[3] == 0;
            case 2:
                return (previous.u[2] == 0) && (previous.u[3] == 0);
            case 1:
            default:
                return (previous.u[1] == 0)
                    && (previous.u[2] == 0)
                    && (previous.u[3] == 0);
            }
        }

        void setParityValue(w128_t parity) {
            state[index] = parity;
            param.parity1 = parity.u[0];
            param.parity2 = parity.u[1];
            param.parity3 = parity.u[2];
            param.parity4 = parity.u[3];
        }

        w128_t getParityValue() const {
            return state[index];
        }

        void setOneBit(int bitPos) {
            setZero();
            int idx = bitPos / 128;
            int p = (bitPos / 64) % 2;
            int r = bitPos % 64;
            state[idx].u[p] = UINT64_C(1) << r;
        }

        /**
         * This method is called by functions in the file
         * simple_shortest_basis.hpp addition of internal state as
         * GF(2) vector is possible when state transition function and
         * output function is GF(2)-linear.
         * @param that SFMT generator added to this generator
         */
        void add(EquidistributionCalculatable<w128_t>& other) {
            sfmt *that = dynamic_cast<sfmt *>(&other);
            if (that == 0) {
                throw std::invalid_argument(
                    "the adder should have same type as the addee.");
            }
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < 2; j++) {
                    state[(i + index) % size].u64[j]
                        ^= that->state[(i + that->index) % size].u64[j];
                }
            }
            previous ^= that->previous;
        }

        int getMexp() const {
            return param.mexp;
        }

        int bitSize() const {
            return size * 128;
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
            reverse_bit_flag = true;
        }

        /**
         * This method is called by the functions in search_temper.hpp
         * to reset the reverse_bit_flag
         */
        void reset_reverse_bit() {
            reverse_bit_flag = false;
        }
        int periodCertification() {
            w128_t tmp;
            w128_t parity;
            parity.u[0] = param.parity1;
            parity.u[1] = param.parity2;
            parity.u[2] = param.parity3;
            parity.u[3] = param.parity4;
            tmp = state[0] & parity;
            int c = count_bit(tmp.u[0]);
            c += count_bit(tmp.u[1]);
            c += count_bit(tmp.u[2]);
            c += count_bit(tmp.u[3]);
            if ((c & 1) == 1) {
                return 1;
            }
            if ((parity.u[0] & 1) == 1) {
                state[0].u[0] ^= 1;
                return 0;
            }
            for (int i = 0; i < 4; i++) {
                uint32_t work = 1;
                for (int j = 0; j < 32; j++) {
                    if ((work & parity.u[i]) != 0) {
                        state[0].u[i] ^= work;
                        return 0;
                    }
                    work = work << 1;
                }
            }
            return 0;
        }
        void d_p() {
            cout << "index = " << dec << index << endl;
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < 4; j++) {
                    cout << hex << setw(8) << setfill('0') << state[i].u[j];
                }
                cout << endl;
            }
        }
    private:
        sfmt& operator=(const sfmt&) {
            throw std::logic_error("can't assign");
        }
        int size;
        int index;
        int start_mode;
        int weight_mode;
        w128_t * state;
        sfmt_param param;
        bool reverse_bit_flag;
        w128_t previous;
    };
}

#endif
