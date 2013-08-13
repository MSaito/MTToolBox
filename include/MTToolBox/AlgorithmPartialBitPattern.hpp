#ifndef MTTOOLBOX_ALGORITHM_PARTIAL_BITPATTERN_HPP
#define MTTOOLBOX_ALGORITHM_PARTIAL_BITPATTERN_HPP
/**
 * @file AlgorithmPartialBitPattern.hpp
 *
 * @brief search tempering parameters so that the random number
 * generator has a good equidistribution property.
 *
 * @author Mutsuo Saito (Hiroshima University)
 * @author Makoto Matsumoto (The University of Tokyo)
 *
 * Copyright (C) 2011 Mutsuo Saito, Makoto Matsumoto,
 * Hiroshima University and The University of Tokyo.
 * All rights reserved.
 *
 * The 3-clause BSD License is applied to this software, see
 * LICENSE.txt
 */
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <unistd.h>
#include <tr1/memory>
#include <MTToolBox/AlgorithmTempering.hpp>
#include <MTToolBox/TemperingCalculatable.hpp>
#include <MTToolBox/AlgorithmEquidistribution.hpp>

namespace MTToolBox {
    /**
     * @class AlgorithmPartialBitPattern
     *
     * @brief 疑似乱数生成器の高次元均等分布性を改善するために、テンパ
     * リングパラメータを探索するアルゴリズム
     *
     * このアルゴリズムはMTGPのテンパリングパラメータを探索するために開
     * 発され、TinyMTのテンパリングパラメータ探索にも使われている。
     *
     * @caution テンパリングパラメータの探索をしても十分良い高次元均等
     * 分布が得られない場合は、状態遷移関数の変更を考慮した方がよい。状
     * 態遷移関数で十分ビットミックスされていない場合、単純なテンパリン
     * グで均等分布次元を最大化することはできないだろう。
     *
     * @tparam T 疑似乱数生成器の出力の型, 例えば uint32_t など。
     *
     * @tparam bit_len テンパリングパラメータのビット長, 通常は出力のビッ
     * ト長と等しいと思われる。
     *
     * @tparam param_num テンパリングパラメータの数, MTGPでは4, TinyMT
     * では 1
     *
     * @tparam try_bit_len 出力のうちテンパリングされる部分の長さ。上位
     * からこのビット数だけがテンパリングされる。
     *
     * @tparam step 一度に何ビットずつビットパターンを生成するか。この
     * 数を大きくすると実行時間が延びるであろう。MTGP では 5 ビットずつ
     * 生成している。
     *
     * @tparam lsb このフラグが指定されると、上位ビットと下位ビットが反
     * 転されて均等分布次元が計算される。つまり下位ビットの均等分布次元
     * をあげたい時に指定する。TestU01のBigCrushには下位ビットの均等分
     * 布次元を改善することによってパス可能性が高まるテストがある。
     */
    template<typename T,
             int bit_len, int param_num, int try_bit_len, int step = 5,
             bool lsb = false>
    class AlgorithmPartialBitPattern : public AlgorithmTempering<T> {
    public:
        /**
         * search tempering parameters.
         * @returns 0
         */
        int operator()(TemperingCalculatable<T>& rand,
                       bool verbose = false) {
            using namespace std;
            if (verbose) {
                cout << "searching..." << endl;
            }
            if (lsb) {
                rand.setReverseOutput();
                if (verbose) {
                    cout << "searching from LSB" << endl;
                }
            } else {
                rand.resetReverseOutput();
                if (verbose) {
                    cout << "searching from MSB" << endl;
                }
            }
            int delta = 1000;
            for (int p = 0; p < try_bit_len; p += step) {
                int max_depth = p + step;
                if (max_depth > try_bit_len) {
                    max_depth = try_bit_len;
                }
                for (int i = 0; i < param_num; i++) {
                    delta = search_best_temper(rand, p, i, max_depth, verbose);
                }
            }
            if (verbose) {
                cout << "delta = " << dec << delta << endl;
            }
            rand.resetReverseOutput();
            return 0;
        }

        bool isLSBTempering() {
            return lsb;
        }
    private:
        void make_temper_bit(TemperingCalculatable<T>& rand,
                             int start, int size,
                             int param_pos, T pattern) {
            T mask = make_mask(start, size);
            rand.setTemperingPattern(mask, pattern, param_pos);
            rand.setUpTempering();
        }

        /**
         * search for one tempering parameter. generate all bit pattern
         * for v-bit to max_v_bit of the tempering parameter, and calculate
         * equidistribution property, then select the best bit pattern.
         * This function calls itself recursively.
         *
         * @param rand linear generator
         * @param v_bit the bit of tempering parameter is searched
         * @param param_pos position of current tempering parameter.
         * @param max_v_bit search stops at this bit.
         * @param verbose if true output searching process
         * @return delta an integer which shows equidistribution property.
         */
        int search_best_temper(TemperingCalculatable<T>& rand, int v_bit,
                               int param_pos, int max_v_bit, bool verbose) {
            using namespace std;
            int bitSize = rand.bitSize();
            int delta;
            int min_delta = bitSize * bit_len;
            T min_pattern = 0;
            int size = max_v_bit - v_bit;
            T pattern;
            T mask = make_mask(v_bit, size);
            //for (int i = 0; i < (1 << size); i++) {
            for (int i = (1 << size) -1; i >= 0; i--) {
                if (lsb) {
                    pattern = static_cast<T>(i) << v_bit;
                } else {
                    pattern = static_cast<T>(i) << (bit_len - v_bit - size);
                }
                rand.setTemperingPattern(mask, pattern, param_pos);
                delta = get_equidist(rand, bit_len, bitSize);
                if (delta < min_delta) {
                    min_delta = delta;
                    min_pattern = pattern;
                } else if (delta == min_delta) {
                    if (count_bit(min_pattern) < count_bit(pattern)) {
                        if (verbose) {
                            cout << "pattern chagne" << hex << min_pattern
                                 << ":" << pattern << endl;
                        }
                        min_pattern = pattern;
                    }
                }
            }
            make_temper_bit(rand, v_bit, size, param_pos, min_pattern);
            if (verbose) {
                cout << min_delta << ":" << hex << min_pattern << endl;
            }
            return min_delta;
        }

        /**
         * wrapper of shortest_basis#get_equidist()
         *
         * @param rand linear generator
         * @param bit_len_ \b bit_len_ MSBs of equidistribution
         * property of the generator is calculated.
         * @param bitSize Mersenne Exponent of the period of rand.
         * @returns summation of equidistribution property
         * from 0 to \b bit_len_ -1 MSBs.
         */
        int get_equidist(TemperingCalculatable<T>& rand,
                         int bit_len_,
                         int bitSize) {
            //TemperingCalculatable<T> r(rand);
            AlgorithmEquidsitribution<T> sb(rand, bit_len_);
            int veq[bit_len_];
            sb.get_all_equidist(veq);
            int sum = 0;
            for (int i = 0; i < bit_len_; i++) {
                sum += (bitSize / (i + 1) - veq[i]) * (bit_len_ - i);
            }
            return sum;
        }
        T get_equidist_pattern(TemperingCalculatable<T>& rand,
                               int bit_len_,
                               int bitSize) {
            //G r(rand);
            T pattern = 0;
            AlgorithmEquidsitribution<T> sb(rand, bit_len_);
            int veq[bit_len_];
            sb.get_all_equidist(veq);
            for (int i = 0; i < bit_len_; i++) {
                if (bitSize / (i + 1) != veq[i]) {
                    pattern |= 1 << (bit_len_ - i - 1);
                }
            }
            return pattern;
        }
        /**
         * make mask which has \b size bits 1s form start.<br/>
         * ex.<br/>
         * when lsb is false, start = 3, size = 5 then
         * return 0x1f000000 (0001 1111 0000 0000 0000 0000 0000 0000)<br/>
         * ex2.<br/>
         * when lsb is true, start = 2, size = 12 then
         * return 0x3ffc (0000 0000 0000 0000 0011 1111 1111 1100)<br/>
         * @param start start bit
         * @param size size of 1s
         * @return bit mask
         */
        T make_mask(int start, int size) {
            T mask = 0;
            mask = ~mask;
            if (start + size > bit_len) {
                size = bit_len - start;
            }
            if (lsb) {
                mask >>= start;
                mask <<= bit_len - size;
                mask >>= bit_len - start - size;
            } else {
                mask <<= start;
                mask >>= bit_len - size;
                mask <<= bit_len - start - size;
            }
            return mask;
        }
    };
}

#endif // MTTOOLBOX_ALGORITHM_PARTIAL_BITPATTERN_HPP
