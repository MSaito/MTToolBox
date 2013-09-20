#ifndef MTTOOLBOX_UTIL_HPP
#define MTTOOLBOX_UTIL_HPP
/**
 * @file util.hpp
 *
 * @brief ユーティリティ関数群
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
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <inttypes.h>
#include <stdint.h>
#include <stdexcept>
#include <tr1/memory>
#include <NTL/GF2X.h>

#if defined(USE_SHA)
#include <openssl/sha.h>
#endif

namespace MTToolBox {
    inline static int count_bit(uint16_t x);
    inline static int count_bit(uint32_t x);
    inline static int count_bit(uint64_t x);
    inline static uint32_t reverse_bit(uint32_t x);
    inline static uint64_t reverse_bit(uint64_t x);

    /**
     * T のビットサイズを返す。
     *
     * sizeof に 8 を掛けているだけなので正確ではない。
     *
     * @tparam T ビットサイズを知りたい型
     * @return T のビットサイズ
     */
    template<typename T>
    int bit_size() {
        return static_cast<int>(sizeof(T) * 8);
    }

    /**
     * 使用しない変数のワーニングを止める
     * @param[in] x 使用しない変数へのポインタ
     */
    inline static void UNUSED_VARIABLE(void * x) {
        (void)x;
    }

    /**
     * n を越えない最大の2のべき乗を返す。
     * @tparam T n の整数型
     * @param[in] n 整数
     * @returns n を越えない最大の2のべき乗
     */
    template<typename T>
    T floor2p(T n) {
        if (n == 1) {
            return 1;
        } else {
            return 2 * floor2p<T>(n / 2);
        }
    }

    /**
     * 出力ストリーム os に多項式 poly を01の文字列で出力する。
     * 次数の低い項の係数が先に出力される。（昇巾順）
     *
     * @param[in,out] os 出力ストリーム
     * @param[in] poly 出力される多項式
     * @param[in] breakline 真なら32文字出力ごとに改行される。
     */
    inline static void print_binary(std::ostream& os,
                                    NTL::GF2X& poly,
                                    bool breakline = true) {
        using namespace NTL;
        if (deg(poly) < 0) {
            os << "0deg=-1" << std::endl;
            return;
        }
        for(int i = 0; i <= deg(poly); i++) {
            if(rep(coeff(poly, i)) == 1) {
                os << '1';
            } else {
                os << '0';
            }
            if (breakline && ((i % 32) == 31)) {
                os << std::endl;
            }
        }
        os << "deg=" << deg(poly) << std::endl;
    }

    /**
     * input を start と end の間の数に変換する。
     * 偏りは気にしない。
     * @param input 入力
     * @param start 範囲の開始
     * @param end 範囲の終了
     * @return \b start <= \b r <= \b end となるような r
     */
    template<typename T>
    int get_range(T input, int start, int end) {
        if (end < start) {
            printf("get_range:%d, %d\n", start, end);
            exit(0);
        }
        return input % (end - start + 1) + start;
    }

    /**
     * GF(2)ベクトルのパラメータテーブルから、より高速で冗長なルックアップテーブルを作成する。
     * @tparam T テーブル内の符号なし整数の型
     * @param dist_tbl 作成されるルックアップテーブル
     * @param src_tbl 元になるGF(2)ベクトルのテーブル
     * @param size \b dist_table　のサイズ
     */
    template<typename T>
    void fill_table(T dist_tbl[], T src_tbl[], int size) {
        for(int i = 1; i < size; i++) {
            for(int j = 1, k = 0; j <= i; j <<= 1, k++) {
                if (i & j) {
                    dist_tbl[i] ^= src_tbl[k];
                }
            }
        }
    }

#if defined(USE_SHA)
    /**
     * GF(2)係数多項式のSHA1ダイジェストを計算する。
     *
     * 多項式の係数を昇巾の順で01からなる文字列に変換し、その後その文字
     * 列のSHA1ダイジェストを計算する。計算されたSHA1ダイジェストは十六進文字列
     * として返される
     *
     * @param[out] str 出力文字列
     * @param[in] poly GF(2)係数多項式
     */
    inline static void poly_sha1(std::string& str, const NTL::GF2X& poly) {
        using namespace NTL;
        using namespace std;
        SHA_CTX ctx;
        SHA1_Init(&ctx);
        if (deg(poly) < 0) {
            SHA1_Update(&ctx, "-1", 2);
        }
        for(int i = 0; i <= deg(poly); i++) {
            if(rep(coeff(poly, i)) == 1) {
                SHA1_Update(&ctx, "1", 1);
            } else {
                SHA1_Update(&ctx, "0", 1);
            }
        }
        unsigned char md[SHA_DIGEST_LENGTH];
        SHA1_Final(md, &ctx);
        stringstream ss;
        for (int i = 0; i < SHA_DIGEST_LENGTH; i++) {
            ss << setfill('0') << setw(2) << hex
               << static_cast<int>(md[i]);
        }
        ss >> str;
    }
#endif

    /**
     * 入力をビット列とみなして最上位の1の位置を0とした最も右側の（下位の）1の位置を返す。
     *
     * 例<ul>
     * <li>calc_1pos(0x8000U) は 0を返す。
     * <li>calc_1pos(0x8002U) は14を返す。
     * <li>calc_1pos(0) は -1 を返す。
     * </ul>
     *
     * このアルゴリズムは以下を参照した
     * @see http://aggregate.org/MAGIC/#Trailing Zero Count
     * @param[in] x 入力
     * @return 最下位の1のあるビットの上位から数えた位置を返す。
     */
    static inline int calc_1pos(uint16_t x)
    {
        if (x == 0) {
            return -1;
        }
        int16_t y = (int16_t)x;
        y = count_bit((uint16_t)((y & -y) - 1));
        return 15 - y;
    }

    /**
     * 入力をビット列とみなして最上位の1の位置を0とした最も右側の（下位の）1の位置を返す。
     *
     * 例<ul>
     * <li>calc_1pos(0x80000000U) は 0を返す。
     * <li>calc_1pos(0x80000002U) は30を返す。
     * <li>calc_1pos(0) は -1 を返す。
     * </ul>
     *
     * このアルゴリズムは以下を参照した
     * @see http://aggregate.org/MAGIC/#Trailing Zero Count
     * @param[in] x 入力
     * @return 最下位の1のあるビットの上位から数えた位置を返す。
     */
    static inline int calc_1pos(uint32_t x)
    {
        if (x == 0) {
            return -1;
        }
        int32_t y = (int32_t)x;
        y = count_bit((uint32_t)(y & -y) - 1);
        return 31 - y;
    }

    /**
     * 入力をビット列とみなして最上位の1の位置を0とした最も右側の（下位の）1の位置を返す。
     *
     * 例<ul>
     * <li>calc_1pos(UINT64_C(0x8000000000000000)) は 0を返す。
     * <li>calc_1pos(UINT64_C(0x8000000000000002)) は62を返す。
     * <li>calc_1pos(0) は -1 を返す。
     * </ul>
     *
     * このアルゴリズムは以下を参照にした
     * @see http://aggregate.org/MAGIC/#Trailing Zero Count
     * @param[in] x 入力
     * @return 最下位の1のあるビットの上位から数えた位置を返す。
     */
    static inline int calc_1pos(uint64_t x)
    {
        if (x == 0) {
            return -1;
        }
        int64_t y = (int64_t)x;
        y = count_bit((uint64_t)(y & -y) - 1);
        return 63 - y;
    }

    /**
     * 1 の個数を数える
     * SIMD within a Register algorithm
     * citing from a website http://aggregate.org/MAGIC/
     *
     * @param x ビットパターン
     * @returns x の中の1の個数
     */
    inline static int count_bit(uint16_t x) {
        x -= (x >> 1) & UINT16_C(0x5555);
        x = ((x >> 2) & UINT16_C(0x3333)) + (x & UINT16_C(0x3333));
        x = ((x >> 4) + x) & UINT16_C(0x0f0f);
        x += (x >> 8);
        return (int)(x & 0x1f);
    }

    /**
     * 1 の個数を数える
     * SIMD within a Register algorithm
     * citing from a website http://aggregate.org/MAGIC/
     */
    /**
     * 1 の個数を数える
     * SIMD within a Register algorithm
     * citing from a website http://aggregate.org/MAGIC/
     *
     * @param x ビットパターン
     * @returns x の中の1の個数
     */
    inline static int count_bit(uint32_t x) {
        x -= (x >> 1) & UINT32_C(0x55555555);
        x = ((x >> 2) & UINT32_C(0x33333333)) + (x & UINT32_C(0x33333333));
        x = ((x >> 4) + x) & UINT32_C(0x0f0f0f0f);
        x += (x >> 8);
        x += (x >> 16);
        return (int)(x & 0x3f);
    }

    /**
     * 1 の個数を数える
     * SIMD within a Register algorithm
     * citing from a website http://aggregate.org/MAGIC/
     *
     * @param x ビットパターン
     * @returns x の中の1の個数
     */
    inline static int count_bit(uint64_t x) {
        x -= (x >> 1) & UINT64_C(0x5555555555555555);
        x = ((x >> 2) & UINT64_C(0x3333333333333333))
            + (x & UINT64_C(0x3333333333333333));
        x = ((x >> 4) + x) & UINT64_C(0x0f0f0f0f0f0f0f0f);
        x += (x >> 8);
        x += (x >> 16);
        x += (x >> 32);
        return (int)(x & 0x7f);
    }

    /**
     * 1 の個数を数える
     * SIMD within a Register algorithm
     * citing from a website http://aggregate.org/MAGIC/
     *
     * @param x ビットパターン
     * @returns x の中の1の個数
     */
    inline static uint32_t reverse_bit(uint32_t x)
    {
        uint32_t y = 0x55555555;
        x = (((x >> 1) & y) | ((x & y) << 1));
        y = 0x33333333;
        x = (((x >> 2) & y) | ((x & y) << 2));
        y = 0x0f0f0f0f;
        x = (((x >> 4) & y) | ((x & y) << 4));
        y = 0x00ff00ff;
        x = (((x >> 8) & y) | ((x & y) << 8));
        return((x >> 16) | (x << 16));
    }

    /**
     * 1 の個数を数える
     * SIMD within a Register algorithm
     * citing from a website http://aggregate.org/MAGIC/
     *
     * @param x ビットパターン
     * @returns x の中の1の個数
     */
    inline static uint64_t reverse_bit(uint64_t x)
    {
        uint64_t y = UINT64_C(0x5555555555555555);
        x = (((x >> 1) & y) | ((x & y) << 1));
        y = UINT64_C(0x3333333333333333);
        x = (((x >> 2) & y) | ((x & y) << 2));
        y = UINT64_C(0x0f0f0f0f0f0f0f0f);
        x = (((x >> 4) & y) | ((x & y) << 4));
        y = UINT64_C(0x00ff00ff00ff00ff);
        x = (((x >> 8) & y) | ((x & y) << 8));
        y = UINT64_C(0x0000ffff0000ffff);
        x = (((x >> 16) & y) | ((x & y) << 16));
        return((x >> 32) | (x << 32));
    }

}
#endif //MTTOOLBOX_UTIL_HPP
