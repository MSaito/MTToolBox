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
        /** 均等分布次元 */
        int kv;
        /** テンパリングパラメータの数 */
        int size;

        /**
         * 単純なコンストラクタ
         * @param[in] param_num テンパリングパラメータの数
         */
        temper_params(int param_num) {
            size = param_num;
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
        /**
         * kv を無視した比較
         * @param that 比較対象
         */
        bool eql(temper_params<T>& that) {
            if (size != that.size) {
                return false;
            }
            bool result = true;
            for (int i = 0; i < size; i++) {
                if (param[i] != that.param[i]) {
                    result = false;
                    break;
                }
            }
            return result;
        }
#if defined(DEBUG)
        string toString() {
            stringstream ss;
            ss << "[";
            ss << dec << kv << ":";
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
         * @tparam T 疑似乱数生成器の出力の型
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
            initial->kv = 0;
            for (int i = 0; i < size; i++) {
                initial->param[i] = 0;
            }

            vector<shared_ptr<tempp> > params;
            params.push_back(initial);
            int kv;
#if defined(DEBUG)
            cout << "DEBUG: bit_len = " << dec << bit_len << endl;
#endif
//            int min_shift = get_min_shift(shifts, size);
            for (int p = 0; p < limit; p++) {
                vector<shared_ptr<tempp> > current;
                current.clear();
                for (unsigned int i = 0; i < params.size(); i++) {
                    search_best_temper(rand, p, *params[i],
                                       current, verbose);
                }
#if defined(DEBUG)
                cout << "DEBUG: current:size = " << dec << current.size()
                     << endl;
                for (unsigned int i = 0; i < current.size(); i++) {
                    cout << "DEBUG:current[" << dec << i << "]:"
                         << current[i]->toString() << endl;
                }
#endif
                kv = -1;
                for (unsigned int i = 0; i < current.size(); i++) {
                    if (current[i]->kv > kv) {
                        kv = current[i]->kv;
                    }
#if defined(DEBUG) && 0
                    else {
                        cout << "DEBUG: current[i]->kv = " << dec
                             << current[i]->kv << endl;
                    }
#endif
                }
#if defined(DEBUG)
                cout << "DEBUG: kv = " << dec << kv << endl;
#endif
                if (verbose) {
                    cout << "kv = " << dec << kv << endl;
                }
                params.clear();
                //T mask = 0;
                //mask = (~mask) << (obSize - p - 1);
                for (unsigned int i = 0; i < current.size(); i++) {
                    if (current[i]->kv == kv) {
                        //params.push_back(current[i]);
#if 0
                        for (int j = 0; j < size; j++) {
                            current[i]->param[j] &= mask;
                        }
#endif
                        push_if_notexist(params, current[i]);
                    }
                }
#if defined(DEBUG)
                cout << "DEBUG: params:size = " << dec << params.size()
                     << endl;
                for (unsigned int i = 0; i < params.size(); i++) {
                    cout << "DEBUG: params[" << dec << i << "]:"
                         << params[i]->toString() << endl;
                }
#endif
            }
#if 0
            // kv を最大にするテンパリングパラメータの中でハミングウェイト最大
            // のものをテンパリングパラメータにする。(複数のなかで最初に見つかったもの)
            int hamming = -1;
            int index = -1;
            for (unsigned int i = 0; i < params.size(); i++) {
                int h = 0;
#if defined(DEBUG)
                cout << "DEBUG:params[" << dec << i << "]"
                     << params[i]->toString() << endl;
#endif
                for (int j = 0; j < size; j++) {
                    h += count_bit(params[i]->param[j]);
                }
                if (hamming < h) {
                    hamming = h;
                    index = i;
                }
            }
#endif
            T mask = 0;
            mask = ~mask;
            for (int i = 0; i < size; i++) {
                //rand.setTemperingPattern(mask, params[index]->param[i], i);
                rand.setTemperingPattern(mask, params[0]->param[i], i);
            }
            rand.setUpTempering();
            if (verbose) {
                cout << "kv = " << dec << kv << endl;
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
//        bool lsb;
        int * shifts;
        int num_pat;
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
        void search_best_temper(TemperingCalculatable<T>& rand,
                                int v_bit,
                                const tempp& para,
                                vector<shared_ptr<tempp> >& current,
                                bool verbose) {
            int kv = -1;
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
                cout << "*pattern:" << pattern->toString() << endl;
#endif
                for (int j = 0; j < size; j++) {
                    rand.setTemperingPattern(mask, pattern->param[j], j);
                }
                pattern->kv = get_equidist(rand, v_bit + 1);
                if (verbose) {
                    cout << "pattern->kv:" << dec << pattern->kv << endl;
                }
                if (pattern->kv >= kv) {
                    push_if_notexist(current, pattern);
                    kv = pattern->kv;
                }
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
        int get_equidist(TemperingCalculatable<T>& rand,
                         int bit_length) {
            AlgorithmEquidistribution<T> sb(rand, bit_length);
            int veq[bit_length];
            sb.get_all_equidist(veq);
            return veq[bit_length - 1];
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
            result.kv = 0;
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
#if defined(DEBUG)
                    cout << "make_pattern:" << dec << pat << ","
                         << dec << v << "," << result.toString()
                         << hex << para_mask
                         << "," << dec << sum << endl;
#endif
                    result.param[index] &= ~(one << (obSize - v - sum - 1));
#if defined(DEBUG)
                    cout << "make_pattern:" << dec << pat << ","
                         << dec << v << "," << result.toString()
                         << hex << para_mask
                         << "," << dec << sum
                         << "," << hex << ~(one << (obSize - v - sum - 1))
                         << "," << dec << v + sum + 1
                         << "," << obSize
                         << endl;
#endif
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
         * @param[in] shifts シフト量
         * @return true 該当する場合
         */
        bool inRange(int pat, int v) const {
#if 0
            for (int i = 0; i < size; i++) {
                if ((obSize < v + shifts[i]) && ((pat >> i) != 0)) {
                    return false;
                }
            }
#endif
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

        int get_min_shift(int nums[], int length) {
            int min = 5000;
            for (int i = 0; i < length; i++) {
                if (min > nums[i]) {
                    min = nums[i];
                }
            }
            return min;
        }

        void push_if_notexist(vector<shared_ptr<tempp> >& vec_temp,
                              shared_ptr<tempp>& pattern) {
            bool found = false;
            for (unsigned int i = 0; i < vec_temp.size(); i++) {
                if (pattern->eql(*vec_temp[i])) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                vec_temp.push_back(pattern);
            }
        }
    };
}

#endif // MTTOOLBOX_ALGORITHM_BEST_BITS_HPP
