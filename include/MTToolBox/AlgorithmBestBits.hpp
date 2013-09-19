#ifndef MTTOOLBOX_ALGORITHM_BEST_BITS_HPP
#define MTTOOLBOX_ALGORITHM_BEST_BITS_HPP
/**
 * @file AlgorithmBestBits.hpp
 *
 * @brief テンパリングパラメータを探索するアルゴリズムであり、MTDCで使
 * 用されているものである。
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
     * @tparam T 個々のテンパリングパラメータの型
     * @tparam size テンパリングパラメータの個数
     */
    template<typename T, int size>
    class temper_params {
    public:
        /** テンパリングパラメータ */
        T param[size];
        /** 均等分布次元の理論値との差の総和 */
        int delta;

        /**
         * 単純なコンストラクタ
         * @param[in] p テンパリングパラメータ
         * @param[in] d 均等分布次元の理論値との差の総和
         */
        temper_params(const T p[], int d) {
            for (int i = 0; i < size; i++) {
                param[i] = p[i];
            }
            delta = d;
        }
    };

    /**
     * @class AlgorithmBestBits
     *
     * @brief テンパリングパラメータを探索するアルゴリズム(MT用)
     *
     * 疑似乱数生成器の高次元均等分布性を改善するために、テンパリングパ
     * ラメータを探索するアルゴリズムこのアルゴリズムはMTDCのテンパリン
     * グパラメータ探索をシミュレートする。
     *
     * \b 注意： テンパリングパラメータの探索をしても十分良い高次元均等
     * 分布が得られない場合は、状態遷移関数の変更を考慮した方がよい。状
     * 態遷移関数で十分ビットミックスされていない場合、単純なテンパリン
     * グで均等分布次元を最大化することはできないだろう。
     *
     * @tparam T 疑似乱数生成器の出力の型, 例えば uint32_t など。
     * @tparam bit_len テンパリングパラメータのビット長, 通常は出力のビット長と等しい
     * と思われる。
     * @tparam param_num テンパリングパラメータの数, MTDCでは2。
     * @tparam shifts テンパリングパラメータと対になるシフト数。正は左シフト。
     * @tparam lsb このアルゴリズムでは使わない
     */
    template<typename T,
             int bit_len,
             int param_num,
             const int shifts[],
             bool lsb = false>
    class AlgorithmBestBits : public AlgorithmTempering<T> {
    public:
        /**
         * テンパリングパラメータのクラス
         */
        typedef temper_params<T, param_num> tempp;

        /**
         * テンパリングパラメータを探索する
         *
         * テンパリングパラメータを生成し、均等分布次元を計算してよいテ
         * ンパリングパラメータを決定する。求められたよいテンパリングパ
         * ラメータは疑似乱数生成器にセットされる。
         *
         * @tparam T 疑似乱数生成器の出力の型
         * @param[in,out] rand 疑似乱数生成器
         * @param[in] verbose 余分な情報を出力するフラグ
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
         * テンパリングパラメータを探索する。
         *
         * @param[in,out] rand GF(2)線形疑似乱数生成器、TemperingSearcherのサブクラス
         * @param[in] v_bit 今からテンパリングしようとするビット(上から数えた)
         * @param[in] para v_bit -1 ビット目までのテンパリングパラメータ
         * @param[out] current v ビット目のテンパリングパラメータとデルタのvector
         * @param[in] shifts シフト量の配列
         * @param[in] verbose true なら余計な情報を標準出力に表示する。
         */
        void search_best_temper(TemperingSearcher<T>& rand,
                                int v_bit,
                                const tempp& para,
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

        /**
         * AlgorithmEquidistribution#get_equidist()のラッパー
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

        /**
         * v ビット目のテンパリングパラメータを作成する。
         *
         * 既に計算済みの MSBから v - 1 ビットのテンパリングパラメータである
         * \b para \b に1ビット分足して MSB から v ビットのテンパリングパラメータを
         * 作成する。
         * @param[out] result 作成されるテンパリングパラメータ
         * @param[in] pat 不明
         * @param[in] v 今回作成するビット
         * @param[in] para これまでに作成したテンパリングパラメータ
         */
        void make_pattern(tempp& result,
                          int pat, int v, const tempp& para) const {
            int obSize = bit_size<T>();
            for (int i = 0; i < size; i++) {
                T mask = 1 & (pat >> i);
                mask = mask << (obSize - 1 - v);
                result.pattern[i] = para.pattern[i] | mask;
            }
        }

        /**
         * マスク作成の該当部分かチェックする。
         *
         * シフトによって0になる部分のビットパターンは作らないので、該
         * 当するかどうかチェックする。
         *
         * @param[in] pat ビットパターン
         * @param[in] v ビット位置
         * @param[in] shifts シフト量
         * @return true 該当する場合
         */
        bool inRange(int pat, int v, const int shifts[]) const {
            int obSize = bit_size<T>();
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
