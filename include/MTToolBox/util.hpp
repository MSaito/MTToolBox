#ifndef MTTOOLBOX_UTIL_HPP
#define MTTOOLBOX_UTIL_HPP
/**
 * @file util.hpp
 *
 *\japanese
 * @brief ユーティリティ関数群
 *\endjapanese
 *
 *\english
 * @brief Utility functions
 *\endenglish
 *
 * @author Mutsuo Saito (Hiroshima University)
 * @author Makoto Matsumoto (Hiroshima University)
 *
 * Copyright (C) 2013, 2016 Mutsuo Saito, Makoto Matsumoto
 * and Hiroshima University.
 * All rights reserved.
 *
 * The 3-clause BSD License is applied to this software, see
 * LICENSE.txt
 */
//#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <inttypes.h>
#include <stdint.h>
#include <stdexcept>
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
     *\japanese
     * T のビットサイズを返す。
     *
     * sizeof に 8 を掛けているだけなので正確ではない。
     *
     * @tparam T ビットサイズを知りたい型
     * @return T のビットサイズ
     *\endjapanese
     *
     *\english
     * Returns bit size of T
     *
     * This is not correct because just multiply sizeof(T) by eight.
     *
     * @tparam T type
     * @return bit size of T
     *\endenglish
     */
    template<typename T>
    int bit_size() {
        return static_cast<int>(sizeof(T) * 8);
    }

    /**
     *\japanese
     * 使用しない変数のワーニングを止める
     * @param[in] x 使用しない変数へのポインタ
     *\endjapanese
     *
     *\english
     * Stop warning of unused variable
     * @param[in] x pointer to unused variable
     *\endenglish
     */
    inline static void UNUSED_VARIABLE(void * x) {
        (void)x;
    }

    /**
     *\japanese
     * n を越えない最大の2のべき乗を返す。
     * @tparam T n の整数型
     * @param[in] n 整数
     * @returns n を越えない最大の2のべき乗
     *\endjapanese
     *
     *\english
     * Return greatest power of two not greater than \b n
     * @tparam T type of \b n
     * @param[in] n integer
     * @return greatest 2<sup>x</sup> \< = n
     *\endenglish
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
     *\japanese
     * 出力ストリーム \b os に多項式 \b poly を01の文字列で出力する。
     *
     * 次数の低い項の係数が先に出力される。（昇巾順）
     *
     * @param[in,out] os 出力ストリーム
     * @param[in] poly 出力される多項式
     * @param[in] breakline 真なら32文字出力ごとに改行される。
     *\endjapanese
     *
     *\english
     * Outputs coefficients of \b poly to \b os.
     *
     * Coefficients of low degree terms are outputted first.
     * @param[in, out] os output stream
     * @param[in] poly polynomial to be outputted
     * @param[in] breakline if true, line break after every 32 terms.
     *\endenglish
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
     *\japanese
     * input を start と end の間の数に変換する。
     *
     * 偏りは気にしない。
     * @param[in] input 入力
     * @param[in] start 範囲の開始
     * @param[in] end 範囲の終了
     * @return \b start \<= \b r \<= \b end となるような r
     *\endjapanese
     *
     *\english
     * Changes input to number between start and end
     *
     * This conversion may not be uniformly distributed.
     * @param[in] input input
     * @param[in] start start of range
     * @param[in] end end of range
     * @return r satisfies that \b start \<= r \<= \b end
     *\endenglish
     */
    template<typename T>
    int get_range(T input, int start, int end) {
        if (end < start) {
            //printf("get_range:%d, %d\n", start, end);
            std::cout << "get_range:" << start << ", " << end << std::endl;
            exit(0);
        }
        return static_cast<int>(input % static_cast<unsigned int>((end - start + 1)))
            + start;
    }

    /**
     *\japanese
     * GF(2)ベクトルのパラメータテーブルから、より高速で冗長なルックアップテーブルを作成する。
     * @tparam T テーブル内の符号なし整数の型
     * @param[out] dist_tbl 作成されるルックアップテーブル
     * @param[in] src_tbl 元になるGF(2)ベクトルのテーブル
     * @param[in] size \b dist_tbl　のサイズ
     *\endjapanese
     *
     *\english
     * Makes a fast and redundant lookup table from a tabale of GF(2) vectors.
     * @tparam T type of elements in tables.
     * @param[out] dist_tbl a lookup table to be made
     * @param[in] src_tbl table of GF(2) vectors.
     * @param[in] size size of \b dist_tbl
     *\endenglish
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
     *\japanese
     * GF(2)係数多項式のSHA1ダイジェストを計算する。
     *
     * 多項式の係数を昇巾の順で01からなる文字列に変換し、その後その文字
     * 列のSHA1ダイジェストを計算する。計算されたSHA1ダイジェストは十六
     * 進文字列として返される
     *
     * @param[out] str 出力文字列
     * @param[in] poly GF(2)係数多項式
     *\endjapanese
     *
     *\english
     * Calculates SHA1 digest of polynomial over GF(2).
     *
     * Convert polynomial to a string of 0 and 1, where coefficient
     * of lowest degree appears first. Then Calculate SHA1 digest of
     * the string. Finally, the calculated digest is converted to
     * a hexadecimal string and returned.
     *
     * @param[out] str output string
     * @param[in] poly polynomial over GF(2)
     *\endenglish
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
     *\japanese
     * 入力をビット列とみなして最上位の1の位置を0とした最も右側の（下位の）1の位置を返す。
     *
     * 1の位置を求めるアルゴリズムは以下のページのものを使用した。
     * @see http://aggregate.org/MAGIC/#Trailing Zero Count
     * @param[in] x 入力
     * @return 最下位の1のあるビットの上位から数えた位置を返す。
     *\endjapanese
     *
     *\english
     * Returns the position of 1 which appears lowest (most right
     * side) in \b x, where the position of MSB becomes zero.
     *
     * Algorithm to search position of 1 is from this page;
     * @see http://aggregate.org/MAGIC/#Trailing Zero Count
     * @param[in] x input
     * @return the position of 1 which appears lowest (most right
     * side) in \b x, where the position of MSB becomes zero.
     *\endenglish
     */
    static inline int calc_1pos(uint16_t x)
    {
        if (x == 0) {
            return -1;
        }
        int16_t y = (int16_t)x;
        y = static_cast<int16_t>(count_bit(static_cast<uint16_t>((y & -y) - 1)));
        return 15 - y;
    }

    /**
     *\japanese
     * \copydoc calc_1pos(uint16_t)
     * 例<ul>
     * <li>calc_1pos(0x80000000U) は 0を返す。
     * <li>calc_1pos(0x80000002U) は30を返す。
     * <li>calc_1pos(0) は -1 を返す。
     * </ul>
     * @param[in] x 入力
     * @return 最下位の1のあるビットの上位から数えた位置を返す。
     *\endjapanese
     *
     *\english
     * \copydoc calc_1pos(uint16_t)
     * Example<ul>
     * <li>calc_1pos(0x80000000U) returns 0.
     * <li>calc_1pos(0x80000002U) returns 30.
     * <li>calc_1pos(0) returns -1.
     * </ul>
     * @param[in] x input
     * @return the position of 1 which appears lowest (most right
     * side) in \b x, where the position of MSB becomes zero.
     *\endenglish
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
     *\japanese
     * \copydoc calc_1pos(uint16_t)
     * 例<ul>
     * <li>calc_1pos(UINT64_C(0x8000000000000000)) は 0を返す。
     * <li>calc_1pos(UINT64_C(0x8000000000000002)) は62を返す。
     * <li>calc_1pos(0) は -1 を返す。
     * </ul>
     * @param[in] x 入力
     * @return 最下位の1のあるビットの上位から数えた位置を返す。
     *\endjapanese
     *
     *\english
     * \copydoc calc_1pos(uint16_t)
     * Example<ul>
     * <li>calc_1pos(UINT64_C(0x8000000000000000)) returns 0.
     * <li>calc_1pos(UINT64_C(0x8000000000000002)) returns 62.
     * <li>calc_1pos(0) returns -1.
     * </ul>
     * @param[in] x input
     * @return the position of 1 which appears lowest (most right
     * side) in \b x, where the position of MSB becomes zero.
     *\endenglish
     */
    static inline int calc_1pos(uint64_t x)
    {
        if (x == 0) {
            return -1;
        }
        int64_t y = (int64_t)x;
        y = count_bit((uint64_t)(y & -y) - 1);
        return 63 - static_cast<int>(y);
    }

    /**
     *\japanese
     * 1 の個数を数える
     *
     * レジスタ内SIMDアルゴリズム
     * http://aggregate.org/MAGIC/ より引用
     *
     * @param[in] x ビットパターン
     * @returns x の中の1の個数
     *\endjapanese
     *
     *\english
     * Counts number of 1s.
     *
     * SIMD within a Register algorithm
     * citing from a website http://aggregate.org/MAGIC/
     * @param[in] x bit pattern
     * @return number of 1s in \b x.
     *\endenglish
     */
    inline static int count_bit(uint16_t x) {
        x -= (x >> 1) & UINT16_C(0x5555);
        x = ((x >> 2) & UINT16_C(0x3333)) + (x & UINT16_C(0x3333));
        x = ((x >> 4) + x) & UINT16_C(0x0f0f);
        x += (x >> 8);
        return (int)(x & 0x1f);
    }

    /**
     *\japanese
     * \copydoc count_bit(uint16_t)
     * @param x ビットパターン
     * @returns x の中の1の個数
     *\endjapanese
     *
     *\english
     * \copydoc count_bit(uint16_t)
     * @param[in] x bit pattern
     * @return number of 1s in \b x.
     *\endenglish
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
     *\japanese
     * \copydoc count_bit(uint16_t)
     * @param x ビットパターン
     * @returns x の中の1の個数
     *\endjapanese
     *
     *\english
     * \copydoc count_bit(uint16_t)
     * @param[in] x bit pattern
     * @return number of 1s in \b x.
     *\endenglish
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
     *\japanese
     * ビットを反転する
     *
     * 入力ビットの上位と下位を反転する。最上位ビットは最下位ビットになる
     * レジスタ内SIMDアルゴリズム
     * http://aggregate.org/MAGIC/ より引用
     * @param x ビットパターン
     * @returns x を反転したビットパターン
     *\endjapanese
     *
     *\english
     * Reverse bits
     *
     * Reverse upper side and lower side in input.
     * The MSB bocomes the LSB.
     * SIMD within a Register algorithm
     * citing from a website http://aggregate.org/MAGIC/
     * @param[in] x bit pattern
     * @return reverse of \b x
     *\endenglish
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
     *\japanese
     * \copydoc reverse_bit(uint32_t)
     * @param x ビットパターン
     * @returns x を反転したビットパターン
     *\endjapanese
     *
     *\english
     * \copydoc count_bit(uint32_t)
     * @param[in] x bit pattern
     * @return reverse of \b x
     *\endenglish
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

    /**
     *\japanese
     * 多項式の最小公倍数を求める。
     * @param[out] lcm x と y の最小公倍数多項式
     * @param[in] x 多項式
     * @param[in] y 多項式
     *\endjapanese
     *
     *\english
     * Calculate the Least Common Multiple of x and y.
     * @param[out] lcm the Least Common Multiple of x and y.
     * @param[in] x a polynomial
     * @param[in] y a polynomial
     *\endenglish
     */
    inline static void LCM(NTL::GF2X& lcm, const NTL::GF2X& x,
                           const NTL::GF2X& y) {
        using namespace NTL;
        GF2X gcd;
        GCD(gcd, x, y);
        mul(lcm, x, y);
        lcm /= gcd;
    }

    /**
     *\japanese
     * 符号なし整数をGF(2)ベクトルに変換する。
     * 上位ビットがベクトルの初めの要素になる。（デバッグの時見やすいように）
     * @tparam U 変換元の型
     * @param[out] result 結果のGF(2)ベクトル
     * @param[in] value 変換元符号なし整数
     *\endjapanese
     *
     *\english
     * convert unsigned integer to GF(2) vector
     * MSB of unsigned integer becomes first element of the vector.
     * @tparam U type of conversion source
     * @param[out] result The result vector
     * @param[in] value source integer
     *\endenglish
     */
    template<typename U>
    inline static void toGF2Vec(NTL::vec_GF2& result, U value) {
        U mask = 1;
        int bitSize = bit_size<U>();
        result.SetLength(bitSize);
        mask = mask << (bitSize - 1);
        for (int i = 0; i < bitSize; i++) {
            if (value & mask) {
                result.put(i, 1);
            } else {
                result.put(i, 0);
            }
            mask = mask >> 1;
        }
    }


    /**
     *\japanese
     * その型の1を返す
     * SIMD型は、そのSIMD型のファイルでこのテンプレートを特殊化する。
     * @tparam U 1を返す型
     * @return 指定された型の1
     *\endjapanese
     *
     *\english
     * return one of specified type
     * @tparam U type of one
     * @return one of specified type
     *\endenglish
     */
    template<typename U>
    inline U getOne() {
        return static_cast<U>(1);
    }

    /**
     *\japanese
     * ゼロをセットする
     * SIMD型は、そのSIMD型のファイルでこのテンプレートを特殊化する。
     * @tparam U 0をセットする型
     * @param[out] x ゼロをセットする変数
     *\endjapanese
     *
     *\english
     * set zero
     * @tparam U type to set zero
     * @param[out] x variable set to be zero
     *\endenglish
     */
    template<typename U>
    inline void setZero(U& x) {
        x = 0;
    }

    /**
     *\japanese
     * 特定位置のビットを求める
     * SIMD型は、そのSIMD型のファイルでこのテンプレートを特殊化する。
     * @tparam U ビットを求める型
     * @param[in] bits ビット列
     * @param[in] pos ビット位置
     * @return 指定位置のビット、ゼロまたは１
     *\endjapanese
     *
     *\english
     * set zero
     * @tparam U type to get bit
     * @param[in] bits variable to get bit
     * @param[in] pos bit position
     * @return bit of position, zero or one.
     *\endenglish
     */
    template<typename U>
    inline unsigned int getBitOfPos(U bits, int pos) {
        return (bits >> pos) & 1;
    }

    /**
     *\japanese
     * 変数の指定位置のビットを1または0にセットする
     * SIMD型は、そのSIMD型のファイルでこのテンプレートを特殊化する。
     * @tparam U 変数の型
     * @param[in,out] bits 指定位置のビットをセットする変数
     * @param[in] pos 指定位置
     * @param[in] b ビット, 0または1
     *\endjapanese
     *
     *\english
     * set zero or 1 to specified bit of specified type variable.
     * @tparam U type of variable
     * @param[in,out] bits variable
     * @param[in] pos position of bit, LSB is 0
     * @param[in] b bit to set
     *\endenglish
     */
    template<typename U>
    inline void setBitOfPos(U *bits, int pos, unsigned int b) {
        b = b & 1;
        U mask = ~(static_cast<U>(1) << pos);
        *bits &= mask;
        *bits |= static_cast<U>(b) << pos;
    }

    /**
     *\japanese
     * GF(2)ベクトルを符号なし整数に変換する。
     * 上位ビットがベクトルの初めの要素になる。（デバッグの時見やすいように）
     * @param[in] value 変換元GF(2)ベクトル
     * @return 変換後符合なし整数
     *\endjapanese
     *
     *\english
     * convert GF(2) vector to unsigned integer
     * the first element of vector becomes the MSB of the result integer.
     * @param[in] value source vector
     * @return resulted integer
     *\endenglish
     */
    template<typename U>
    inline static U fromGF2Vec(NTL::vec_GF2& value) {
        U result;
        setZero(result);
        U mask;
        setZero(mask);
        int bitSize = bit_size<U>();
        int pos = bitSize - 1;
        setBitOfPos(&mask, pos, 1);
        for (int i = 0; i < bitSize; i++) {
            if (!IsZero(value[i])) {
                result |= mask;
            }
            pos--;
            setZero(mask);
            setBitOfPos(&mask, pos, 1);
        }
        return result;
    }

    /**
     *\japanese
     * ゼロかどうか、判定する。
     * SIMD型は、そのSIMD型のファイルでこのテンプレートを特殊化する。
     * @tparam U 変数の型
     * @param[in] x 変数
     * @return 変数がゼロかどうか ゼロならtrue
     *\endjapanese
     *
     *\english
     * check if varible is zero or not
     * @tparam U type of variable
     * @param[in] x variable to be checked
     * @return ture if x is zero
     *\endenglish
     */
    template<typename U>
    inline bool isZero(U x) {
        return x == 0;
    }

    /**
     *\japanese
     * V 型をU型に変換する
     * SIMD型は、そのSIMD型のファイルでこのテンプレートを特殊化する。
     * @tparam U 変換先の型
     * @tparam V 変換元変数の型
     * @param[in] x 変数
     * @return U 型の値
     *\endjapanese
     *
     *\english
     * convert to type U from type V
     * @tparam U type convert to
     * @tparam V type convert from
     * @param[in] x variable
     * @return value of type U
     *\endenglish
     */
    template<typename U, typename V>
    inline U convert(V x) {
        return static_cast<U>(x);
    }
}
#endif //MTTOOLBOX_UTIL_HPP
