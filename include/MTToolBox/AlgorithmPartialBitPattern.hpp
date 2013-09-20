#ifndef MTTOOLBOX_ALGORITHM_PARTIAL_BITPATTERN_HPP
#define MTTOOLBOX_ALGORITHM_PARTIAL_BITPATTERN_HPP
/**
 * @file AlgorithmPartialBitPattern.hpp
 *
 * @brief 疑似乱数生成器の高次元均等分布性を改善するために、テンパ
 * リングパラメータを探索するアルゴリズム
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
     * \b 注意: テンパリングパラメータの探索をしても十分良い高次元均等
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
         * テンパリングパラメータを探索する。
         *
         * この処理の呼び出しは長時間かかる可能性がある。
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

        /**
         * このテンパリングパラメータ探索がLSBからの探索であるかどうか
         *@return true LSBからのテンパリングパラメータ探索
         */
        bool isLSBTempering() const {
            return lsb;
        }
    private:
        void make_temper_bit(TemperingCalculatable<T>& rand,
                             T mask,
                             int param_pos,
                             T pattern) {
            rand.setTemperingPattern(mask, pattern, param_pos);
            rand.setUpTempering();
        }

        /**
         * テンパリングパラメータをひとつ探索する。
         *
         * # v_bit から max_v_bit までのすべてのビットパターンを生成し
         * # すべてのビットパターンについてvビット精度均等分布次元を計算し
         * # Δの最も小さかったビットパターンを選択する。
         * # 同じΔであればハミングウェイトの大きいビットパターンを選択する
         * なお、このメソッドは再帰する。
         *
         * @param rand GF(2)疑似乱数生成器
         * @param v_bit テンパリングパラメータを探索するビット範囲
         * @param param_pos 何番目のテンパリングパラメータを探索しているか
         * @param max_v_bit このビットで探索をやめる
         * @param verbose true なら探索過程の情報を出力する
         * @return delta このテンパリングパラメータの設定結果によるΔ
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
            int length = bit_size<T>() / 4;
            for (int i = (1 << size) -1; i >= 0; i--) {
                if (lsb) {
                    pattern = static_cast<T>(i) << v_bit;
                } else {
                    pattern = static_cast<T>(i) << (bit_len - v_bit - size);
                }
                make_temper_bit(rand, mask, param_pos, pattern);
                delta = get_equidist(rand, bit_len);
                if (delta < min_delta) {
                    if (verbose) {
                        cout << "pattern chagne " << hex << min_pattern
                             << ":" << pattern << endl;
                    }
                    min_delta = delta;
                    min_pattern = pattern;
                } else if (delta == min_delta) {
                    if (count_bit(min_pattern) < count_bit(pattern)) {
                        if (verbose) {
                            cout << "pattern chagne " << hex << min_pattern
                                 << ":" << pattern << endl;
                        }
                        min_pattern = pattern;
                    }
                }
            }
            make_temper_bit(rand, mask, param_pos, min_pattern);
            if (verbose) {
                cout << dec << min_delta << ":"
                     << hex << setw(length) << min_pattern << ":"
                     << hex << setw(length) << mask << endl;
            }
            return min_delta;
        }

        /**
         * AlgorithmEquidistributionの呼び出し
         *
         * @param rand GF(2)疑似乱数生成器
         * @param bit_len_ MSB から \b bit_len_ 分の均等分布次元を計算する
         * @returns Δ 理論値と実現値の差の合計
         */
        int get_equidist(TemperingCalculatable<T>& rand,
                         int bit_len_) {
            AlgorithmEquidistribution<T> sb(rand, bit_len_);
            int veq[bit_len_];
            int sum = sb.get_all_equidist(veq);
            return sum;
        }

        /**
         * ビットマスク作成
         *
         * これはテンパリングパラメータではなく、テンパリングパラメータを
         * セットする部分をマスクするためのビットマスクである。
         *
         * @param start ビットマスクを開始するビット位置
         * @param size ビットマスクの長さ
         * @return ビットマスク
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
