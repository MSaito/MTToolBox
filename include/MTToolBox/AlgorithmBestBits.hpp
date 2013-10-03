#ifndef MTTOOLBOX_ALGORITHM_BEST_BITS_HPP
#define MTTOOLBOX_ALGORITHM_BEST_BITS_HPP
/**
 * @file AlgorithmBestBits.hpp
 *
 * @brief テンパリングパラメータを探索するアルゴリズムであり、MTDCで使
 * 用されているものである。
 *
 * @author Mutsuo Saito (Hiroshima University)
 * @author Makoto Matsumoto (Hiroshima University)
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
#include <vector>
#include <MTToolBox/AlgorithmTempering.hpp>
#include <MTToolBox/TemperingCalculatable.hpp>
#include <MTToolBox/AlgorithmEquidistribution.hpp>
#include <MTToolBox/util.hpp>

#if defined(DEBUG)
#include <sstream>
#endif

namespace MTToolBox {
    using namespace std;
    using namespace std::tr1;

    /**
     * @class temper_params
     * @brief テンパリングパラメータのクラス
     *
     * @tparam T 個々のテンパリングパラメータの型
     */
    template<typename T>
    class temper_params {
    public:
        /** テンパリングパラメータ */
        T * param;
        /** 均等分布次元の理論値との差の総和 */
        int delta;
        /** テンパリングパラメータの数 */
        int size;

        /**
         * 単純なコンストラクタ
         * @param[in] param_num テンパリングパラメータの数
         */
        temper_params(int param_num) {
            size = param_num;
            delta = 0;
            param = new T[param_num];
            for (int i = 0; i < param_num; i++) {
                param[i] = 0;
            }
        }

        /**
         * デストラクタ
         */
        ~temper_params() {
            delete[] param;
        }
#if defined(DEBUG)
        string toString() {
            stringstream ss;
            ss << "[";
            ss << dec << delta << ":";
            for (int i = 0; i < size; i++) {
                ss << hex << setw(8) << setfill('0') << param[i] << ",";
            }
            ss << "]";
            return ss.str();
        }
#endif
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
     */
    template<typename T>
    class AlgorithmBestBits : public AlgorithmTempering<T> {
    public:
        /**
         * テンパリングパラメータのクラス
         */
        typedef temper_params<T> tempp;

        /**
         * @param[in] out_bit_length テンパリングパラメータのビット長,
         * 通常は出力のビット長と等しいと思われる。
         * @param[in] shift_values テンパリングパラメータと対になるシフト数。
         * 正は左シフト。現状では右シフトには対応していない。
         * @param[in] param_num テンパリングパラメータの数, MTDCでは2。
         * 7を越えないこと。実際には 2 の場合しかテストしていない。
         * @param[in] limit_v テンパリングパラメータを上位何ビットまでテンパリングするか。
         * ただし、この値+以後のシフト量だけテンパリングする。
         */
        AlgorithmBestBits(int out_bit_length,
                          const int shift_values[],
                          int param_num,
                          int limit_v) {
            limit = limit_v;
            obSize = bit_size<T>();
            bit_len = out_bit_length;
            size = param_num;
            shifts = new int[size];
            num_pat = size * (size + 1) / 2;
            for (int i = 0; i < size; i++) {
                shifts[i] = shift_values[i];
            }
        }

        /**
         * デストラクタ
         */
        ~AlgorithmBestBits() {
            delete[] shifts;
        }

        /**
         * テンパリングパラメータを探索する
         *
         * テンパリングパラメータを生成し、均等分布次元を計算してよいテ
         * ンパリングパラメータを決定する。求められたよいテンパリングパ
         * ラメータは疑似乱数生成器にセットされる。
         *
         * @param[in,out] rand 疑似乱数生成器
         * @param[in] verbose 余分な情報を出力するフラグ
         * @returns 0
         */
        int operator()(TemperingCalculatable<T>& rand,
                       bool verbose = false) {
            rand.resetReverseOutput();
            if (verbose) {
                cout << "searching from MSB" << endl;
            }
            shared_ptr<tempp> initial(new tempp(size));
            initial->delta = 0;
            for (int i = 0; i < size; i++) {
                initial->param[i] = 0;
            }

            vector<shared_ptr<tempp> > params;
            params.push_back(initial);
            int delta;
#if defined(DEBUG)
            cout << "DEBUG: bit_len = " << dec << bit_len << endl;
#endif
            for (int p = 0; p < limit; p++) {
                vector<shared_ptr<tempp> > current;
                current.clear();
                for (unsigned int i = 0; i < params.size(); i++) {
                    search_best_temper(rand, p, *params[i],
                                       current, verbose);
                }
                delta = rand.bitSize() * obSize;
                for (unsigned int i = 0; i < current.size(); i++) {
                    if (current[i]->delta < delta) {
                        delta = current[i]->delta;
                    }
                }
#if defined(DEBUG)
                cout << "DEBUG: delta = " << dec << delta << endl;
#endif
                if (verbose) {
                    cout << "delta = " << dec << delta << endl;
                }
                params.clear();
                for (unsigned int i = 0; i < current.size(); i++) {
                    if (current[i]->delta == delta) {
                        params.push_back(current[i]);
                    }
                }
            }
            T mask = 0;
            mask = ~mask;
            for (int i = 0; i < size; i++) {
                //rand.setTemperingPattern(mask, params[index]->param[i], i);
                rand.setTemperingPattern(mask, params[0]->param[i], i);
            }
            rand.setUpTempering();
            if (verbose) {
                cout << "delta = " << dec << delta << endl;
            }
            rand.resetReverseOutput();
            return 0;
        }

        bool isLSBTempering() {
            return false;
        }
    private:
        int limit;
        int bit_len;
        int size;
        int obSize;
        int * shifts;
        int num_pat;

        /**
         * テンパリングパラメータを探索する。
         *
         * @param[in,out] rand GF(2)線形疑似乱数生成器、TemperingSearcherのサブクラス
         * @param[in] v_bit 今からテンパリングしようとするビット(上から数えた)
         * @param[in] para v_bit -1 ビット目までのテンパリングパラメータ
         * @param[out] current v ビット目のテンパリングパラメータとデルタのvector
         * @param[in] verbose true なら余計な情報を標準出力に表示する。
         */
        void search_best_temper(TemperingCalculatable<T>& rand,
                                int v_bit,
                                const tempp& para,
                                vector<shared_ptr<tempp> >& current,
                                bool verbose) {
            int delta = rand.bitSize() * obSize;
            T mask = 0;
            mask = ~mask;
            // size が 2 なら 11, 10, 01, 00 の4パターン
            // というのが大きな間違いで、実は 8 パターン
            num_pat = size * (size + 1) / 2;
            for (int32_t i = (1 << num_pat) -1; i >= 0; i--) {
                if (! inRange(i, v_bit)) {
                    continue;
                }
                shared_ptr<tempp> pattern(new tempp(size));
                make_pattern(*pattern, i, v_bit, para);
#if defined(DEBUG)
                cout << "pattern:" << pattern->toString() << endl;
#endif
                for (int j = 0; j < size; j++) {
                    rand.setTemperingPattern(mask, pattern->param[j], j);
                }
                pattern->delta = get_equidist(rand, v_bit + 1);
                if (verbose) {
                    cout << "pattern->delta:" << dec << pattern->delta << endl;
                }
                if (pattern->delta <= delta) {
                    current.push_back(pattern);
                    delta = pattern->delta;
                }
            }
        }

        /**
         * AlgorithmEquidistribution#get_equidist()のラッパー
         *
         * @param rand 疑似乱数生成器
         * @param bit_len_ \b bit_len_ 長の MSBs の均等分布次元を計算する
         * @returns 均等分布次元の理論値との差を0からbit_lenまで合計したもの
         */
        int get_equidist(TemperingCalculatable<T>& rand,
                         int bit_length) {
            AlgorithmEquidistribution<T> sb(rand, bit_length);
            int veq[bit_length];
            return sb.get_all_equidist(veq);
        }

        /**
         * v ビット目のテンパリングパラメータを作成する。
         *
         * 既に計算済みの MSBから v - 1 ビットのテンパリングパラメータである
         * \b para \b に1ビット分足して MSB から v ビットのテンパリングパラメータを
         * 作成する。
         * @param[out] result 作成されるテンパリングパラメータ
         * @param[in] pat テンパリングパラメータ数からなるビットパターン
         * @param[in] v 今回作成するビット
         * @param[in] para これまでに作成したテンパリングパラメータ
         */
        void make_pattern(tempp& result,
                          int pat, int v, const tempp& para) const {
            result.delta = 0;
            for (int i = 0; i < result.size; i++) {
                result.param[i] = para.param[i];
            }
            T para_mask = 0;
            para_mask = (~para_mask) >> v;
            int index = 0;
            int idx = 0;
            int rdx = size - 1;
            T mask = 1 << (num_pat - 1);
            int sum = 0;
            T one = 1;
            while (mask != 0) {
                if ((pat & mask) && (obSize > v + sum + 1)) {
                    result.param[index] |= (one << (obSize - v - sum - 1))
                    & para_mask;
                } else if (((pat & mask) == 0) && (obSize > v + sum + 1)) {
                    result.param[index] &= ~(one << (obSize - v - sum - 1));
                }
                mask = mask >> 1;
                index = index + 1;
                if (index >= size) {
                    sum += shifts[rdx];
                    index = idx;
                    idx++;
                    rdx--;
                }
            }

#if defined(DEBUG)
            cout << "make_pattern:" << dec << pat << ","
                 << dec << v << "," << result.toString()
                 << hex << para_mask << endl;
#endif
        }

        /**
         * マスク作成の該当部分かチェックする。
         *
         * シフトによって0になる部分のビットパターンは作らないので、該
         * 当するかどうかチェックする。
         *
         * @param[in] pat ビットパターン
         * @param[in] v ビット位置
         * @return true 該当する場合
         */
        bool inRange(int pat, int v) const {
            int index = 0;
            int idx = 0;
            int rdx = size - 1;
            T mask = 1 << (num_pat - 1);
            int sum = 0;
            while (mask != 0) {
                if ((pat & mask) && (v + shifts[index] + sum > obSize)) {
                    return false;
                }
                mask = mask >> 1;
                index = index + 1;
                if (index >= size) {
                    sum += shifts[rdx];
                    index = idx;
                    idx++;
                    rdx--;
                }
            }
            return true;
        }
    };
}

#endif // MTTOOLBOX_ALGORITHM_BEST_BITS_HPP
