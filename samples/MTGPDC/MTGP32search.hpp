#ifndef MTGP32SEARCH_HPP
#define MTGP32SEARCH_HPP
/**
 * @file MTGP32search.hpp
 *
 * @brief Mersenne Twister for Graphic Processors 32 bit
 * 32 bit class for parameter search
 *
 * @author Mutsuo Saito (Hiroshima University)
 * @author Makoto Matsumoto (Hiroshima University)
 *
 * Copyright (c) 2010 Mutsuo Saito, Makoto Matsumoto and Hiroshima
 * University. All rights reserved.
 * Copyright (c) 2011 Mutsuo Saito, Makoto Matsumoto, Hiroshima
 * University and University of Tokyo. All rights reserved.
 *
 * The new BSD License is applied to this software, see LICENSE.txt
 */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cerrno>
#include <cstring>
#include <sstream>
//#include <unistd.h>
#include <MTToolBox/EquidistributionCalculatable.hpp>
#include <MTToolBox/AlgorithmPartialBitPattern.hpp>
#include <MTToolBox/ParameterGenerator.hpp>
#include <MTToolBox/MersenneTwister.hpp>
#include <MTToolBox/util.hpp>
#include "mtgp_param.h"

/**
 * @namespace mtgp
 * @brief namespace for MTGP
 */
namespace mtgp {
    using namespace MTToolBox;
    using namespace NTL;
    using namespace std;

    /**
     * @class mtgp32
     * @brief Mersenne Twister for Graphic Processors 32 bit
     * class for parameter search. Class mtgp does not have
     * internal state. linear_generator class has internal
     * state and handling functions.
     */
    class mtgp32 : public TemperingCalculatable<uint32_t> {
    public:
        /**
         * constructor with mexp and id
         * @param mexp_ Mersenne exponent
         * @param id_ parameter to generate diffrent sequence
         */
        mtgp32(int mexp_, uint32_t id_) {
            param.mexp = mexp_;
            unsigned long ss = static_cast<unsigned long>(mexp_) /
                (sizeof(uint32_t) * 8) + 1UL;
            state_size = static_cast<int>(ss);
            param.mask = (~UINT32_C(0))
                << (bit_size<uint32_t>() * state_size - param.mexp);
            param.id = id_;
            state = new uint32_t[static_cast<unsigned long>(state_size)];
            idx = 0;
            reverse_bit_flag = false;
            seed(0);
        }

        /**
         * copy constructor
         * @param src source object
         */
        mtgp32(const mtgp32& src) : TemperingCalculatable<uint32_t>(),
                                    param(src.param) {
            state_size = src.state_size;
            state = new uint32_t[static_cast<unsigned long>(state_size)];
            idx = src.idx;
            reverse_bit_flag = src.reverse_bit_flag;
            for (int i = 0; i < state_size; i++) {
                state[i] = src.state[i];
            }
        }

        ~mtgp32() {
            delete[] state;
        }

        mtgp32 * clone() const {
            return new mtgp32(*this);
        }

        /**
         * accessor to Mersenne exponet
         * @returns Mersenne exponet
         */
        int bitSize() const {
            return param.mexp;
        }

        uint32_t generate() {
            next_state();
            return temper();
        }

        /**
         * This method is called by the functions in simple_shortest_basis.hpp
         * This method returns \b bit_len bits of MSB of generated numbers
         * If reverse_bit_flag is set, bits are taken from LSB
         * @param bit_len bit length from MSB or LSB
         * @return generated numbers of bit_len
         */
        uint32_t generate(int bit_len) {
            uint32_t u;
            if (reverse_bit_flag) {
                u = reverse_bit(generate());
            } else {
                u = generate();
            }
            uint32_t mask = 0;
            mask = (~mask) << (32 - bit_len);
            return u & mask;
        }

#if 0
        /**
         * generate block of numbers,
         * designed for fast calculation of minimum polynomial.
         * index of state array is not used because it is not needed.
         * Do not call me twice before parameter set up.
         * duplicate call will return the same output.
         * @param[out] array outputs are placed here
         * @param[in] size number of output
         * @param state internal state array
         * @param state_size size of state
         */
        void block_generate_without_tempering(uint32_t array[],
                                              int size,
                                              uint32_t state[],
                                              int state_size) {
            int i;
            for (i = 0; i < state_size - param.pos; i++) {
                array[i] = rec(state[i], state[i + 1], state[i + param.pos]);
            }
            for (; i < state_size - 1; i++) {
                array[i] = rec(state[i],
                               state[i + 1],
                               array[i + param.pos - state_size]);
            }
            array[state_size - 1] = rec(state[state_size - 1],
                                         array[0],
                                         array[param.pos - 1]);
            for (i = 0; i < size - state_size; i++) {
                array[i + state_size] = rec(array[i],
                                             array[i + 1],
                                             array[i + param.pos]);
            }
        }
#endif
        /**
         * set up parameters randomly.
         * @param mt random number generator
         * @param state_size size of internal state array
         */
        void setUpParam(ParameterGenerator& mt){
           if (state_size > 100) {
                param.pos
                    = get_range(mt.getUint32(), 3,
                                state_size - floor2p<int>(state_size) - 1);
            } else {
                param.pos = get_range(mt.getUint32(), 3, state_size / 2 - 1);
            }
            //param.sh1 = get_range(mt.next(), 1, 31);
            //param.sh2 = get_range(mt.next(), 1, 31);
            param.sh1 = 13;
            param.sh2 = 4;
            for (int i = 0; i < 4; i++) {
                param.tbl[i] = mt.getUint32();
                if ((param.tbl[i] & 0x0f) == (UINT32_C(1) << i)) {
                    param.tbl[i] ^= (UINT32_C(1) << i);
                }
            }
            param.tbl[2] = (param.tbl[2] & UINT32_C(0xfff0000f))
                ^ ((param.id << 4) & UINT32_C(0x000ffff0));
            param.tbl[3] = (param.tbl[3] & UINT32_C(0x0000ffff))
                ^ (param.id & UINT32_C(0xffff0000));
            memset(param.p, 0, sizeof(param.p));
            fill_table(param.p, param.tbl, 16);
        }
        /**
         * returns parameter
         * @return parameter of me
         */
        mtgp_param<uint32_t> get_param() {
            return param;
        }
        const std::string getHeaderString() {
            return param.getHeaderString();
        }

        const std::string getParamString() {
            return param.getParamString();
        }
        void setTemperingPattern(uint32_t mask, uint32_t pattern, int src_bit) {
            param.tmp_tbl[src_bit] &= ~mask;
            param.tmp_tbl[src_bit] |= pattern & mask;
        }
        /**
         * This method is called by functions in the file
         * simple_shortest_basis.hpp addition of internal state as
         * GF(2) vector is possible when state transition function and
         * output function is GF(2)-linear.
         * @param that tinymt generator added to this generator
         */
        void add(EquidistributionCalculatable<uint32_t>& other) {
            mtgp32* that = dynamic_cast<mtgp32 *>(&other);
            if (that == 0) {
                throw std::invalid_argument(
                    "the adder should have the same type as the addee.");
            }
#if defined(DEBUG)
            if (!param.equals(that->param)) {
                cerr << "in add: the adder should have the same parameter"
                     << " as the addee." << endl;
                cerr << "this:";
                param.printParam(cerr);
                cerr << endl;
                cerr << "that:";
                that->param.printParam(cerr);
                throw std::invalid_argument(
                    "the adder should have the same parameter as the addee.");
            }
#endif
            for (int i = 0; i < state_size; i++) {
                state[(i + idx) % state_size]
                    ^= that->state[(i + that->idx) % state_size];
            }
        }

        void setUpTempering() {
            memset(param.tp, 0, sizeof(param.tp));
            fill_table(param.tp, param.tmp_tbl, 16);
        }
        void set_param(mtgp_param<uint32_t> src) {
            param = src;
        }

        /**
         * initialize the internal state array
         * @param seed seed of initialization
         */
        void seed(uint32_t value) {
            state[0] = value;
            for (int i = 1; i < state_size; i++) {
                state[i] = UINT32_C(1812433253)
                    * (state[i - 1] ^ (state[i - 1] >> 30))
                    + static_cast<uint32_t>(i);
            }
            idx = state_size - 1;
        }

        /**
         * This method is called by the functions in simple_shortest_basis.hpp
         */
        void setZero() {
            for (int i = 0; i < state_size; i++) {
                state[i] = 0;
            }
        }

        /**
         * This method is called by the functions in the file
         * simple_shortest_basis.hpp
         * @return true if all elements of status is zero
         */
        bool isZero() const {
            if ((state[idx] & mask) != 0) {
                return false;
            }
            for (int i = 1; i < state_size; i++) {
                if (state[(idx + i) % state_size] != 0) {
                    return false;
                }
            }
            return true;
        }
        /**
         * This method is called by the functions in search_temper.hpp
         * to calculate the equidistribution properties from LSB
         */
        void setReverseOutput() {
            reverse_bit_flag = true;
        }

        /**
         * This method is called by the functions in search_temper.hpp
         * to reset the reverse_bit_flag
         */
        void resetReverseOutput() {
            reverse_bit_flag = false;
        }

        bool isReverseOutput() {
            return reverse_bit_flag;
        }
    private:
        uint32_t mask;
        int state_size;
        mtgp_param<uint32_t> param;
        int idx;
        bool reverse_bit_flag;
        uint32_t * state;

        /**
         * transform internal state
         * @param state internal state array
         * @param state_size size of state
         * @param index of state array
         */
        void next_state() {
            idx = (idx + 1) % state_size;
            state[idx] = rec(state[idx],
                             state[(idx + 1) % state_size],
                             state[(idx + param.pos) % state_size]);
        }
        /**
         */
        uint32_t temper() {
            uint32_t v = state[idx];
            uint32_t t = state[(idx + param.pos - 1) % state_size];
            t = t ^ (t >> 16);
            t = t ^ (t >> 8);
            return v ^ param.tp[t & 0x0f];
        }
        /**
         * make look up tabel of F2 vectors
         * @param dist_tbl look up table to be made, size of \b dist_tbl
         * should be 2<sup>\b size</sup>.
         * @param src_tbl source table
         * @param size element number of \b src_tbl
         */
        void fill_table(uint32_t dist_tbl[], uint32_t src_tbl[], int size) {
            dist_tbl[0] = 0;
            for(int i = 1; i < size; i++) {
                dist_tbl[i] = 0;
                for(int j = 1, k = 0; j <= i; j <<= 1, k++) {
                    if (i & j) {
                        dist_tbl[i] ^= src_tbl[k];
                    }
                }
            }
        }

        /**
         * calculate recursion formula
         * @param p1 further point
         * @param p2 next further point
         * @param p3 middle pickup piont
         * @return new value
         */
        uint32_t rec(uint32_t p1, uint32_t p2, uint32_t p3) {
            p1 = (p1 & param.mask) ^ p2;
            p1 ^= p1 << param.sh1;
            p3 = p1 ^ (p3 >> param.sh2);
            return p3 ^ param.p[p3 & 0x0f];
        }


    };
}

#endif //MTGP32SEARCH_HPP
