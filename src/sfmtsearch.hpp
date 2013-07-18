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
 * @author Mutsuo Saito (Hiroshima University)
 * @author Makoto Matsumoto (The University of Tokyo)
 *
 * Copyright (C) 2013 Mutsuo Saito, Makoto Matsumoto
 * and Hiroshima University.
 * All rights reserved.
 *
 * The 3-clause BSD License is applied to this software, see
 * LICENSE.txt
 */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cerrno>
#include <sstream>
#include <unistd.h>
#include <search_temper.hpp>
#include "mt19937ar.h"

/**
 * @namespace sfmt
 * name space for sfmt
 */
namespace sfmt {
    using namespace MTToolBox;
    using namespace NTL;
    using namespace std;

    union w128_t {
        vector unsigned int s;
        uint32_t u[4];
        uint64_t u64[2];
    };

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

        /**
         * This method is used in output.hpp.
         * @return header line of output.
         */
        string get_header() {
            return "mexp, pos1, sl1, sl2, sr1, sr2, msk1, msk2, msk3, msk4";
        }

        /**
         * This method is used in output.hpp.
         * @return string of parameters
         */
        string get_string() {
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
            string s;
            ss >> s;
            return s;
        }

        /**
         * This method is used for DEBUG.
         * @return string of parameters.
         */
        string get_debug_string() {
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
    class sfmt {
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
            reverse_bit_flag = false;
        }

        ~sfmt() {
            delete[] state;
        }

        /**
         * The copy constructor.
         * @param src The origin of copy.
         */
        sfmt(const sfmt& src) : param(src.param) {
            size = src.param.mexp / 128 + 1;
            state = new w128_t[size];
            reverse_bit_flag = false;
        }

        /**
         * Constructor by parameter.
         * @param src_param
         */
        sfmt(const sfmt_param& src_param) : param(src_param) {
            size = src_param.mexp / 128 + 1;
            state = new w128_t[size];
            reverse_bit_flag = false;
        }

        /**
         * This method is called by the functions in simple_shortest_basis.hpp
         */
        void make_zero_state() {
            for (int i = 0; i < size; i++) {
                state[i] = 0;
            }
            index = 0;
        }

        /**
         * This method returns mexp
         * @return mexp
         */
        int get_mexp() const {
            return param.mexp;
        }

        /**
         * This method returns size of internal state
         * @return size
         */
        int get_state_size() const {
            return size;
        }

        /**
         * This method initialize internal state.
         * This initialization is simple.
         * @param seed seed for initialization
         */
        void seeding(uint32_t seed) {
            make_zero_state();
            if (param.pos1 > 0 && param.pos1 < size) {
                state[param.pos1 - 1] = seed;
            } else {
                state[0] = seed;
            }
            index = 0;
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
            do_recursion(&state[index],
                         &state[index],
                         &state[(index + param.pos1) % size],
                         &state[(index + size - 2) % size],
                         &state[(index + size -2) % size]);
            index = (index + 1) % size;
        }

        /**
         * get a part of internal state without tempering
         * @return
         */
        w128_t get_uint() {
            return state[index];
        }

        /**
         * getter of state
         * @param index index of internal state
         * @return element of internal state at \b index
         */
        w128_t get_state(int index) {
            return state[index];
        };

        /**
         * Important method, generate new random number
         * @return new pseudo random number
         */
        w128_t generate() {
            next_state();
            return state[index];
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
        void setup_param(uint32_t num) {
            uint32_t work = num ^ (num << 15) ^ (num << 23);
            work <<= 1;
            param.mat1 = (work & 0xffff0000) | (param.id & 0xffff);
            param.mat2 = (work & 0xffff) | (param.id & 0xffff0000);
            param.mat1 ^= param.mat1 >> 19;
            param.mat2 ^= param.mat2 << 18;
            param.mat2 ^= 1;
            param.tmat[0] = 0;
        };

        /**
         * getter of parameter
         * @return parameter
         */
        const sfmt_param& get_param() const {
            return param;
        };

        /**
         * getter of recursion parameter
         * @param mat1 parameter mat1 is set after calling this method
         * @param mat2 parameter mat2 is set after calling this method
         */
        void get_mat(uint32_t *mat1, uint32_t *mat2) const {
            *mat1 = param.mat1;
            *mat2 = param.mat2;
        };

        /**
         * This method gives information to the functions in the file
         * search_temper.hpp
         * @return always 1
         */
        int get_temper_param_num() const {
            return 1;
        };

        /**
         * This method is called by the functions in the file
         * simple_shortest_basis.hpp
         * @return true if all elements of state is zero
         */
        bool is_zero() {
            return (state[0] == 0) &&
                (state[1] == 0) &&
                (state[2] == 0) &&
                (state[3] == 0);
        };

        /**
         * This is important,
         * This function tries to improve output quality of randomness.
         * One important point of SFMT is this tempering function,
         * which is not GF(2)-linear.
         * But in calculating parameter phase, NON_LINEAR is never defined.
         * @return improved random number
         */
        uint32_t temper() {
#if defined(NO_TEMPER)
            return state[0];
#else
            uint32_t t0, t1;
#if defined(NON_LINEAR)
            t0 = state[3];
            t1 = state[0] + (state[2] >> sh8);
#else
            t0 = state[3];
            t1 = state[0] ^ (state[2] >> sh8);
#endif
            t0 ^= t1;
            if (t1 & 1) {
                t0 ^= param.tmat[0];
            }
            return t0;
#endif
        };

        /**
         * This method is called by functions in the file search_temper.hpp
         * @param mask available bits of pattern
         * @param pattern bit pattern
         * @param src_bit only 0 is allowed
         */
        void set_temper_pattern(uint32_t mask, uint32_t pattern, int src_bit) {
            param.tmat[src_bit] &= ~mask;
            param.tmat[src_bit] |= pattern & mask;
        };

        /**
         * This method is called by functions in the file
         * simple_shortest_basis.hpp addition of internal state as
         * GF(2) vector is possible when state transition function and
         * output function is GF(2)-linear.
         * @param that SFMT generator added to this generator
         */
        void add(sfmt& that) {
            state[0] ^= that.state[0];
            state[1] ^= that.state[1];
            state[2] ^= that.state[2];
            state[3] ^= that.state[3];
        };

        /**
         * This method is called by functions in the file search_temper.hpp
         * Do not remove this.
         */
        void setup_temper() {
        };

        /**
         * setter of parameter
         * @param src new parameter
         */
        void set_param(sfmt_param src) {
            param = src;
        };

        /**
         * output parameters
         * @param out output stream
         */
        void out_param(ostream& out) {
            string s = param.get_debug_string();
            out << s << endl;
        };

        /**
         * This method is called by the functions in search_temper.hpp
         * to calculate the equidistribution properties from LSB
         */
        void set_reverse_bit() {
            reverse_bit_flag = true;
        };

        /**
         * This method is called by the functions in search_temper.hpp
         * to reset the reverse_bit_flag
         */
        void reset_reverse_bit() {
            reverse_bit_flag = false;
        };
    private:
        int size;
        int index;
        uint32_t * state;
        sfmt_param param;
        bool reverse_bit_flag;
    };
}

#endif
