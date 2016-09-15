#ifndef MTTOOLBOX_ALGORITHM_BEST_BITS_HPP
#define MTTOOLBOX_ALGORITHM_BEST_BITS_HPP
/**
 * @file AlgorithmBestBits.hpp
 *
 *\japanese
 * @brief テンパリングパラメータを探索するアルゴリズムであり、MTDC [1]で使
 * 用されているものである。
 *\endjapanese
 *\english
 *
 * @brief Algorithm for searching tempering parameters.
 * Dynamic Creator[1] uses this algorithm.
 *
 *\endenglish
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
 *
 * Reference:
 * [1] M. Matsumoto and T. Nishimura,
 * Dynamic Creation of Pseudorandom number generator,
 * Monte Carlo and Quasi-Monte Carlo Methods 1998,
 * Springer-Verlag, 56-69, 2000
 * http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/DC/dc.html
 */
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <unistd.h>
#if __cplusplus >= 201103L
#include <memory>
#else // memory
#define MTTOOLBOX_USE_TR1
#include <tr1/memory>
#endif
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
#if defined(MTTOOLBOX_USE_TR1)
    using std::tr1::shared_ptr;
#else
    using std::shared_ptr;
#endif

    /**
     * @class temper_params
     *\japanese
     * @brief テンパリングパラメータのクラス
     *
     * @tparam U 個々のテンパリングパラメータの型, 符号なし型でなければならない
     *\endjapanese
     *\english
     * @brief class which keeps tempering parameters
     *
     * @tparam U type of tempering parameters, should be unsigned type.
     *\endenglish
     */
    template<typename U>
    class temper_params {
    public:
        /**
         *\japanese
         * テンパリングパラメータの配列
         *\endjapanese
         *\english
         * array of tempering parameters
         *\endenglish
         */
        U * param;
        /**
         *\japanese
         * 均等分布次元の理論値との差の総和
         *\endjapanese
         *\english
         * Sum of differences between theoretical upper bound and
         * realized value of dimension of equi-distribution.
         *\endenglish
         */
        int delta;
        /**
         *\japanese
         * テンパリングパラメータの数
         *\endjapanese
         *\english
         * Number of tempering parameters.
         *\endenglish
         */
        int size;

        /**
         *\japanese
         * 単純なコンストラクタ
         * @param[in] param_num テンパリングパラメータの数
         *\endjapanese
         *\english
         * Simple constructor.
         * @param[in] param_num number of tempering parameters
         *\endenglish
         */
        temper_params(int param_num) {
            size = param_num;
            delta = 0;
            param = new U[param_num];
            for (int i = 0; i < param_num; i++) {
                param[i] = 0;
            }
        }

        /**
         *\japanese
         * デストラクタ
         *\endjapanese
         *\english
         * Simple destructor
         *\endenglish
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
     *\japanese
     * @brief テンパリングパラメータを探索するアルゴリズム(MT用)
     *
     * 疑似乱数生成器の高次元均等分布性を改善するために、テンパリングパ
     * ラメータを探索するアルゴリズム。このアルゴリズムはMTDCのテンパリン
     * グパラメータ探索をシミュレートする。
     *
     * \warning テンパリングパラメータの探索をしても十分良い高次元均等
     * 分布が得られない場合は、状態遷移関数の変更を考慮した方がよい。状
     * 態遷移関数で十分ビットミックスされていない場合、単純なテンパリン
     * グで均等分布次元を最大化することはできないだろう。
     *
     * @tparam U 疑似乱数生成器の出力の型, 符号なし型であること、例えば uint32_t など。
     *\endjapanese
     *\english
     * @brief Algorithm which searches tempering parameters.
     *
     * An algorithm which searches tempering parameters which
     * improve equi-distribution property of pseudo random number
     * generator. This algorithm simulate the algorithm of
     * Dynamic creator of Mersenne Twister.
     *
     * \warning If you could not get high dimension of
     * equi-distribution after tempering using parameters got by this
     * algorithm, you should consider changing the design of random
     * number generation algorithm.
     *
     * @tparam U type of output of pseudo random number generator.
     * Should be unsigned type.
     *\endenglish
     */
    template<typename U, typename V = U>
    class AlgorithmBestBits : public AlgorithmTempering<U, V> {
    public:
        /**
         *\japanese
         * テンパリングパラメータのクラス
         *\endjapanese
         *\english
         * a class which keeps tempering parameters.
         *\endenglish
         */
        typedef temper_params<U> tempp;

        /**
         *\japanese
         * @param[in] out_bit_length テンパリングパラメータのビット長,
         * 通常は出力のビット長と等しいと思われる。
         * @param[in] shift_values テンパリングパラメータと対になるシフト数。
         * 正は左シフト。現状では右シフトには対応していない。
         * @param[in] param_num テンパリングパラメータの数, MTDCでは2。
         * 7を越えないこと。実際には 2 の場合しかテストしていない。
         * この値を大きくすると使用メモリおよび実行時間が著しく増大するであろう。
         * @param[in] limit_v テンパリングパラメータを上位何ビットまでテンパリングするか。
         * ただし、この値+以後のシフト量だけテンパリングする。この値を大きくすると
         * 使用メモリおよび実行時間が著しく増大するであろう。
         *\endjapanese
         *\english
         * @param[in] out_bit_length Bit length of tempering parameters.
         * usually, this equals to bit length of an output.
         * @param[in] shift_values Shift values corresponds to
         * tempering parameter. Positive integers mean left shift and
         * negative integers are not supported in current version.

         * @param[in] param_num Number of tempering parameters.  In
         * MTDC, this is two. Do not specify numbers over 7.  Test is
         * only done for 2. Greater numbers will consume huge CPU time
         * and huge memories.

         * @param[in] limit_v limit of tempering bit. Each tempering
         * parameters are searched for limit bit plus shift values of
         * parameters following to the parameters. Greater limit will
         * consume huge CPU time and huge memories.
         *\endenglish
         */
        AlgorithmBestBits(int out_bit_length,
                          const int shift_values[],
                          int param_num,
                          int limit_v) {
            limit = limit_v;
            obSize = bit_size<U>();
            bit_len = out_bit_length;
            size = param_num;
            shifts = new int[size];
            num_pat = size * (size + 1) / 2;
            for (int i = 0; i < size; i++) {
                shifts[i] = shift_values[i];
            }
        }

        /**
         *\japanese
         * デストラクタ
         *\endjapanese
         *\english
         * A destructor
         *\endenglish
         */
        ~AlgorithmBestBits() {
            delete[] shifts;
        }

        /**
         *\japanese
         * テンパリングパラメータを探索する
         *
         * テンパリングパラメータを生成し、均等分布次元を計算してよいテ
         * ンパリングパラメータを決定する。求められたよいテンパリングパ
         * ラメータは疑似乱数生成器にセットされる。
         *
         * @param[in,out] rand 疑似乱数生成器
         * @param[in] verbose 余分な情報を出力するフラグ
         * @returns 0
         *\endjapanese
         *\english
         *
         * Search tempering parameters.
         *
         * Generate tempering parameters and Calculate dimension of
         * equi-distribution and then choose tempering parameters
         * which gives high dimension of equi-distribution.  Searched
         * tempering parameters are set to the pseudo random number
         * generator.
         *
         * @param[in,out] rand A pseudo random number generator
         * @param[in] verbose A flag for printing redundant message
         * @return 0 always zero
         *\endenglish
         */
        int operator()(TemperingCalculatable<U, V>& rand,
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
            int delta = 0;
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
            U mask = 0;
            mask = ~mask;
            for (int i = 0; i < size; i++) {
                rand.setTemperingPattern(mask, params[0]->param[i], i);
            }
            rand.setUpTempering();
            if (verbose) {
                cout << "delta = " << dec << delta << endl;
            }
            rand.resetReverseOutput();
            return 0;
        }

        /**
         *\japanese
         * LSBからのテンパリングかを示す。常にfalse。
         * @return 常にfalse
         *\endjapanese
         *\english
         * Shows if tempering is from LSB.
         * @return Always false.
         *\endenglish
         */
        bool isLSBTempering() const {
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
         *\japanese
         * テンパリングパラメータを探索する。
         *
         * @param[in,out] rand GF(2)線形疑似乱数生成器、TemperingSearcherのサブクラス
         * @param[in] v_bit 今からテンパリングしようとするビット(上から数えた)
         * @param[in] para v_bit -1 ビット目までのテンパリングパラメータ
         * @param[out] current v ビット目のテンパリングパラメータとデルタのvector
         * @param[in] verbose true なら余計な情報を標準出力に表示する。
         *\endjapanese
         *\english
         *\endenglish
         */
        void search_best_temper(TemperingCalculatable<U, V>& rand,
                                int v_bit,
                                const tempp& para,
                                vector<shared_ptr<tempp> >& current,
                                bool verbose) {
            int delta = rand.bitSize() * obSize;
            U mask = 0;
            mask = ~mask;
            // size が 2 なら 111, 110, 101, 100, 011, 010, 001, 000 の8パターン
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
         *\japanese
         * AlgorithmEquidistribution#get_equidist()のラッパー
         *
         * @param rand 疑似乱数生成器
         * @param bit_len_ \b bit_len_ 長の MSBs の均等分布次元を計算する
         * @returns 均等分布次元の理論値との差を0からbit_lenまで合計したもの
         *\endjapanese
         *\english
         * A wrapper of AlgorithmEquidistribution#get_equidist().
         *
         * @param rand A pseudo random number generator.
         * @param bit_length Calculate dimension of equi-distribution
         * of \b bit_length MSBs.
         * @returns A sum of differences of dimension of equi-distribution
         * between theoretical upper bound and realized value.
         *\endenglish
         */
        int get_equidist(TemperingCalculatable<U, V>& rand,
                         int bit_length) {
            AlgorithmEquidistribution<U, V> sb(rand, bit_length);
            int veq[bit_length];
            return sb.get_all_equidist(veq);
        }

        /**
         *\japanese
         * v ビット目のテンパリングパラメータを作成する。
         *
         * 既に計算済みの MSBから v - 1 ビットのテンパリングパラメータである
         * \b para \b に1ビット分足して MSB から v ビットのテンパリングパラメータを
         * 作成する。
         * @param[out] result 作成されるテンパリングパラメータ
         * @param[in] pat テンパリングパラメータ数からなるビットパターン
         * @param[in] v 今回作成するビット
         * @param[in] para これまでに作成したテンパリングパラメータ
         *\endjapanese
         *\english
         *
         * Make a part of tempering parameters \b v-th bit from MSB.
         *
         * \b para keeps tempering parameters which have \b v - 1 bits
         * from MSBs. This method make \b v-th bits of tempering
         * parameters and (\b v + shift_value)-th bits.
         *
         * @param[out] result Created tempering parameters.
         * @param[in] pat Bits pattern used to generate tempering parameters.
         * @param[in] v A bit position to be created.
         * @param[in] para One of previously generated tempering parameters.
         *
         *\endenglish
         */
        void make_pattern(tempp& result,
                          int pat, int v, const tempp& para) const {
            result.delta = 0;
            for (int i = 0; i < result.size; i++) {
                result.param[i] = para.param[i];
            }
            U para_mask = 0;
            para_mask = (~para_mask) >> v;
            int index = 0;
            int idx = 0;
            int rdx = size - 1;
            U mask = 1 << (num_pat - 1);
            int sum = 0;
            U one = 1;
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
         *\japanese
         * マスク作成の該当部分かチェックする。
         *
         * シフトによって0になる部分のビットパターンは作らないので、該
         * 当するかどうかチェックする。
         *
         * @param[in] pat ビットパターン
         * @param[in] v ビット位置
         * @return true 該当する場合
         *\endjapanese
         *\english
         *
         * Check if pattern should be processed.
         *
         * Bit pattern which becomes zero by shift out need not
         * to be created.
         *
         * @param[in] pat Bit pattern
         * @param[in] v Bit position
         * @return true if bit should be created.
         *
         *\endenglish
         */
        bool inRange(int pat, int v) const {
            int index = 0;
            int idx = 0;
            int rdx = size - 1;
            U mask = 1 << (num_pat - 1);
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
#if defined(MTTOOLBOX_USE_TR1)
#undef MTTOOLBOX_USE_TR1
#endif
#endif // MTTOOLBOX_ALGORITHM_BEST_BITS_HPP
