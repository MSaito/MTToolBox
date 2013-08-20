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
    class mtgp32 {
    public:
        /**
         * constructor with mexp and id
         * @param mexp_ Mersenne exponent
         * @param id_ parameter to generate diffrent sequence
         */
        mtgp32(int mexp_, uint32_t id_) {
            memset(&param, 0, sizeof(param));
            param.mexp = mexp_;
            status_size = mexp_ / (sizeof(uint32_t) * 8) + 1;
            param.mask = (~static_cast<uint32_t>(0))
                << (sizeof(uint32_t) * 8 * status_size - param.mexp);
            param.id = id_;
        }
        /**
         * copy constructor
         * @param src source object
         */
        mtgp32(const mtgp32& src) : param(src.param) {
            status_size = src.status_size;
        }
        /**
         * accessor to Mersenne exponet
         * @returns Mersenne exponet
         */
        int get_mexp() {
            return param.mexp;
        }
        /**
         * accessor to status_size
         * @returns size of internal state.
         */
        int get_status_size() {
            return status_size;
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
        /**
         * transform internal state
         * @param status internal state array
         * @param status_size size of status
         * @param index of state array
         */
        void next_state(uint32_t status[], int status_size,
                        int idx) {
            status[idx] = rec(status[idx],
                              status[(idx + 1) % status_size],
                              status[(idx + param.pos) % status_size]);
        }
        /**
         * generate block of numbers,
         * designed for fast calculation of minimum polynomial.
         * index of state array is not used because it is not needed.
         * Do not call me twice before parameter set up.
         * duplicate call will return the same output.
         * @param[out] array outputs are placed here
         * @param[in] size number of output
         * @param status internal state array
         * @param status_size size of status
         */
        void block_generate_without_tempering(uint32_t array[],
                                              int size,
                                              uint32_t status[],
                                              int status_size) {
            int i;
            for (i = 0; i < status_size - param.pos; i++) {
                array[i] = rec(status[i], status[i + 1], status[i + param.pos]);
            }
            for (; i < status_size - 1; i++) {
                array[i] = rec(status[i],
                               status[i + 1],
                               array[i + param.pos - status_size]);
            }
            array[status_size - 1] = rec(status[status_size - 1],
                                         array[0],
                                         array[param.pos - 1]);
            for (i = 0; i < size - status_size; i++) {
                array[i + status_size] = rec(array[i],
                                             array[i + 1],
                                             array[i + param.pos]);
            }
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
         * set up parameters randomly.
         * @param mt random number generator
         * @param status_size size of internal state array
         */
        void setup_param(AbstractGenerator<uint32_t>& mt, int status_size){
            if (status_size > 100) {
                param.pos
                    = get_range(mt.next(), 3,
                                status_size - floor2p<int>(status_size) - 1);
            } else {
                param.pos = get_range(mt.next(), 3, status_size / 2 - 1);
            }
            //param.sh1 = get_range(mt.next(), 1, 31);
            //param.sh2 = get_range(mt.next(), 1, 31);
            param.sh1 = 13;
            param.sh2 = 4;
            for (int i = 0; i < 4; i++) {
                param.tbl[i] = mt.next();
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
        /**
         * returns number of tempering parameter vectors
         * @return number of tempering parameter vectors
         */
        int get_temper_param_num() {
            return 4;
        }
        /**
         */
        uint32_t temper(uint32_t v, uint32_t status[],
                        int status_size, int idx) {
            uint32_t t = status[(idx + param.pos - 1) % status_size];
            t = t ^ (t >> 16);
            t = t ^ (t >> 8);
            return v ^ param.tp[t & 0x0f];
        }
        void set_temper_bit(int target_bit, int src_bit) {
            param.tmp_tbl[src_bit] |= static_cast<uint32_t>(1) << target_bit;
        }
        void reset_temper_bit(int target_bit, int src_bit) {
            param.tmp_tbl[src_bit] &= ~(static_cast<uint32_t>(1) << target_bit);
        }
        void set_temper_pattern(uint32_t mask, uint32_t pattern, int src_bit) {
            param.tmp_tbl[src_bit] &= ~mask;
            param.tmp_tbl[src_bit] |= pattern & mask;
        }
        void setup_temper() {
            memset(param.tp, 0, sizeof(param.tp));
            fill_table(param.tp, param.tmp_tbl, 16);
        }
        void set_param(mtgp_param<uint32_t> src) {
            param = src;
        }
        void out_param(ostream& out) {
            out << "pos:" << param.pos << endl;
            out << "sh1:" << param.sh1 << endl;
            out << "sh2:" << param.sh2 << endl;
            for (int i = 0; i < 4; i++) {
                out << "tbl[" << dec << i << "]:0x" << hex
                    << setw(sizeof(uint32_t) * 2) << setfill('0')
                    << param.tbl[i]
                    << endl;
            }
            for (int i = 0; i < 4; i++) {
                out << "tmp_tbl[" << dec << i << "]:0x"
                    << hex << setw(sizeof(uint32_t) * 2) << setfill('0')
                    << param.tmp_tbl[i] << endl;
            }
            out << "id:" << "0x" << hex << setw(sizeof(uint32_t) * 2)
                << setfill('0') << param.id << dec << endl;
        }
    private:
        int status_size;
        mtgp_param<uint32_t> param;
    };
}

#endif //MTGP32SEARCH_HPP
