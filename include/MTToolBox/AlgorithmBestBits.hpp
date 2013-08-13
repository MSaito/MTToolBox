#ifndef MTTOOLBOX_ALGORITHM_BEST_BITS_HPP
#define MTTOOLBOX_ALGORITHM_BEST_BITS_HPP
/**
 * @file AlgorithmBestBits.hpp
 *
 * @brief テンパリングパラメータを探索するアルゴリズムであり、MTDCで使用されているものである。
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
#include <unistd.h>
#include <tr1/memory>
#include <MTToolBox/AlgorithmTempering.hpp>
#include <MTToolBox/TemperingCalculatable.hpp>
#include <MTToolBox/AlgorithmEquidistribution.hpp>

namespace MTToolBox {
    /**
     * @class temper_params
     * @brief テンパリングパラメータのクラス
     *
     */
    template<typename T, int size>
    class temper_params {
    public:
        T param[size];
        int delta;
        temper_params(T& p, int d) {
            for (int i = 0; i < size; i++) {
                param[i] = p[i];
            }
            delta = d;
        }
    };

    /**
     * @class AlgorithmPartialBitPattern
     *
     * @brief 疑似乱数生成器の高次元均等分布性を改善するために、テンパリングパラメータを
     * 探索するアルゴリズム
     *
     * このアルゴリズムはMTDCのテンパリングパラメータ探索をシミュレートする。
     *
     * @caution テンパリングパラメータの探索をしても十分良い高次元均等分布が得られない
     * 場合は、状態遷移関数の変更を考慮した方がよい。状態遷移関数で十分ビットミックスされて
     * いない場合、単純なテンパリングで均等分布次元を最大化することはできないだろう。
     *
     * @tparam T 疑似乱数生成器の出力の型, 例えば uint32_t など。
     * @tparam bit_len テンパリングパラメータのビット長, 通常は出力のビット長と等しい
     * と思われる。
     * @tparam param_num テンパリングパラメータの数, MTDCでは2。
     * @tparam shifts テンパリングパラメータと対になるシフト数。正は左シフト。
     * @tparam lsb このフラグが指定されると、上位ビットと下位ビットが反転されて均等分布
     * 次元が計算される。つまり下位ビットの均等分布次元をあげたい時に指定する。
     * TestU01のBigCrushには下位ビットの均等分布次元を改善することに
     * よってパスする可能性が高まるテストがある。
     */
    template<typename T,
             int bit_len,
             int param_num,
             const int shifts[],
             bool lsb = false>
    class AlgorithmBestBits : public AlgorithmTempering<T> {
    public:
        typedef temper_params<T, param_num> tempp;

        /**
         * search tempering parameters.
         * @returns 0
         */
        int operator()(TemperingSearcher<T>& rand,
                       bool verbose = false) {
            using namespace std;
            if (verbose) {
                cout << "searching..." << endl;
            }
            if (lsb) {
                rand.setReverseBit();
                if (verbose) {
                    cout << "searching from LSB" << endl;
                }
            } else {
                rand.resetReverseBit();
                if (verbose) {
                    cout << "searching from MSB" << endl;
                }
            }
            shared_ptr<tempp> initial(new tempp());
            initial->delta = 0;
            for (int i = 0; i < param_num; i++) {
                initial->param[i] = 0;
            }

            vector<shared_ptr<tempp> > params();
            params.push_back(initial);
            int delta;
            for (int p = 0; p < bit_len; p++) {
                vector<shared_ptr<tempp> > current();
                current.clear();
                for (int i = 0; i < params.size(); i++) {
                    search_best_temper(rand, p, params[i],
                                       current, shifts, verbose);
                }
                delta = rand.bitSize();
                for (int i = 0; i < current.size(); i++) {
                    if (current[i]->delta < delta) {
                        delta = current[i]->delta;
                    }
                }
                if (verbose) {
                    cout << "delta = " << dec << delta << endl;
                }
                params.clear();
                for (int i = 0; i < current.size(); i++) {
                    if (current[i]->delta == delta) {
                        params.push_back(current[i]);
                    }
                }
            }
            int hamming = 0;
            int index = -1;
            for (int i = 0; i < params.size(); i++) {
                int h = 0;
                for (int j = 0; j < size; j++) {
                    h += bit_count(params[i]->param[j]);
                }
                if (hamming < h) {
                    hamming = h;
                    index = i;
                }
            }
            T mask = 0;
            mask = ~mask;
            for (int i = 0; i < size; i++) {
                rand.setTemperingPattern(mask, params[index]->param[i], i);
            }
            rand.setUpTempering();
            if (verbose) {
                cout << "delta = " << dec << delta << endl;
            }
            rand.resetReverseBit();
            return 0;
        }
        bool isLSBTempering() {
            return lsb;
        }
    private:
        /**
         * search for one tempering parameter. generate all bit pattern
         * for v-bit to max_v_bit of the tempering parameter, and calculate
         * equidistribution property, then select the best bit pattern.
         * This function calls itself recursively.
         *
         * @param rand GF(2)線形疑似乱数生成器、TemperingSearcherのサブクラス
         * @param v_bit 今からテンパリングしようとするビット(上から数えた)
         * @param tempp v_bit -1 ビット目までのテンパリングパラメータ
         * @param current v ビット目のテンパリングパラメータとデルタのvector
         * @param verbose true なら余計な情報を標準出力に表示する。
         */
        void search_best_temper(TemperingSearcher<T>& rand,
                                int v_bit,
                                tempp& para,
                                vector<shared_ptr<tempp> >& current,
                                const int shifts[],
                                bool verbose) {
            using namespace std;
            int bitSize = rand.bitSize();
            int delta;
            T mask = 0;
            mask = ~mask;
            for (int i = (1 << size) -1; i >= 0; i--) {
                if (! inRange(i, v_bit, shifts)) {
                    continue;
                }
                shared_ptr<tempp> pattern(new tempp);
                make_pattern(*pattern, i, v_bit, para);
                for (int j = 0; j < size; j++) {
                    rand.setTemperingPattern(mask, pattern->pattern[j], j);
                }
                delta = get_equidist(rand, v_bit, bitSize);
                pattern->delta = delta;
                current.push_back(pattern);
            }
        }
#if 0
        void make_temper_bit(TemperingSearcher<T>& rand, int start, int size,
                             int param_pos, T pattern) {
            T mask = make_mask(start, size);
            rand.setTemperingPattern(mask, pattern, param_pos);
            rand.setUpTempering();
        }
#endif
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
        int get_equidist(TemperingSearcher<T>& rand,
                         int bit_len_,
                         int bitSize) {
            TemperingSearcher<T> r(rand);
            calc_equidist<TemperingSearcher<T>, T> sb(r, bit_len_);
            int veq[bit_len_];
            sb.get_all_equidist(veq);
            int sum = 0;
            for (int i = 0; i < bit_len_; i++) {
                sum += (bitSize / (i + 1) - veq[i]) * (bit_len_ - i);
            }
            return sum;
        }
#if 0
        T get_equidist_pattern(TemperingSearcher<T>& rand,
                               int bit_len_,
                               int bitSize) {
            G r(rand);
            T pattern = 0;
            shortest_basis<G, T> sb(r, bit_len_);
            int veq[bit_len_];
            sb.get_all_equidist(veq);
            for (int i = 0; i < bit_len_; i++) {
                if (bitSize / (i + 1) != veq[i]) {
                    pattern |= 1 << (bit_len_ - i - 1);
                }
            }
            return pattern;
        }
#endif
        /**
         * 既に計算済みの MSBから v - 1 ビットのテンパリングパラメータである
         * \code para に1ビット分足して MSB から v ビットのテンパリングパラメータを
         * 作成する。
         * @param[out] result 作成されるテンパリングパラメータ
         * @param size size of 1s
         * @return bit mask
         */
        void make_pattern(tempp& result, int pat, int v, tempp& para) const {
            int obSize = bit_size(T);
            for (int i = 0; i < size; i++) {
                T mask = 1 & (pat >> i);
                mask = mask << (obSize - 1 - v);
                result.pattern[i] = para.pattern[i] | mask;
            }
        }

        /**
         * シフトによって0になる部分のビットパターンは作らないので、
         * 該当するかどうかチェックする。
         *
         * \example
         */
        bool inRange(int pat, int v, const int shifts[]) const {
            int obSize = bit_size(T);
            for (int i = 0; i < size; i++) {
                if ((obSize < v + shifts[i]) && ((pat >> i) != 0)) {
                    return false;
                }
            }
            return true;
        }
    };
}

#endif // MTTOOLBOX_ALGORITHM_BEST_BITS_HPP
