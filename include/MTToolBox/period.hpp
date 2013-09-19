#ifndef MTTOOLBOX_PERIOD_HPP
#define MTTOOLBOX_PERIOD_HPP
/**
 * @file period.hpp
 *
 * @brief 最小多項式の計算と原始性の判定
 *
 * @author Mutsuo Saito (Hiroshima University)
 * @author Makoto Matsumoto (The University of Tokyo)
 *
 * Copyright (C) 2013 Mutsuo Saito, Makoto Matsumoto,
 * Hiroshima University and The University of Tokyo.
 * All rights reserved.
 *
 * The 3-clause BSD License is applied to this software, see
 * LICENSE.txt
 */

#include <stdint.h>
#include <NTL/GF2X.h>
#include <NTL/vector.h>
#include <MTToolBox/AbstractGenerator.hpp>

namespace MTToolBox {
    /**
     * 最小多項式を求める
     *
     * @tparam U 疑似乱数生成器の出力の型
     * @param[out] poly 最小多項式
     * @param[in] generator GF(2)疑似乱数生成器
     * @param[in] pos 出力の下位から何ビット目を見るかを指定する
     */
    template<typename U> void
    minpoly(NTL::GF2X& poly, AbstractGenerator<U>& generator, int pos = 0)
    {
        using namespace std;
        using namespace NTL;

        Vec<GF2> v;
        int size = generator.bitSize();
        v.SetLength(2 * size);
        for (int i = 0; i < 2 * size; i++) {
            v[i] = (generator.generate() >> pos) & 1;
        }
        MinPolySeq(poly, v, size);
    }

    /**
     * 2^degree -1 が素数となるかどうかを返す
     *
     * メルセンヌ素数の指数のリストを元に判定しているので完全ではない。
     * 疑似乱数生成器として使用する範囲であればカバーしている。
     * @param[in] degree 判定するべき数
     * @return true 2<sup>degree</sup> -1 が素数の場合
     */
    bool isMexp(uint32_t degree);

    /**
     * 既約判定
     *
     * @param[in] poly GF(2)係数多項式
     * @return true poly が既約の場合
     */
    bool isIrreducible(const NTL::GF2X& poly);

    /**
     * 原始性判定
     *
     * この原始性判定は簡易版であり、poly の次数がメルセンヌ指数の場合のみ
     * 正しい結果を返す。状態空間のビットサイズがメルセンヌ指数でない場合は、
     * この関数を使うべきではない。
     *
     * @param[in] poly GF(2)係数多項式
     * @return true poly の次数がメルセンヌ指数で、かつpolyが既約のとき
     */
    bool isPrime(const NTL::GF2X& poly);

    /**
     * 原始性判定
     *
     * この関数は、poly の次数が degree でない場合は false を返す。
     * prime_factors には2<sup>degree</sup>-1 の素因数分解に現れる素数の
     * リストを与える。2<sup>3</sup>のように通常は同じ素数を複数含むが、
     * 多重度は考慮せずにひとつの素数を１回だけ含むリストを与えればよい。
     * prime_factors が正しくないと結果も正しくないであろう。
     *
     * @param[in] poly GF(2)係数多項式
     * @param[in] degree polyに期待する次数
     * @param[in] prime_factors 2<sup>degree</sup>-1 の素因数分解に現れる素数のリスト
     * @return true poly が原始多項式の場合
     */
    bool isPrime(const NTL::GF2X& poly, int degree,
                 const NTL::Vec<NTL::ZZ>& prime_factors);

    /**
     * 原始性判定
     *
     * この関数は、poly の次数が degree でない場合は false を返す。
     * prime_factors には2<sup>degree</sup>-1 の素因数分解に現れる素数の
     * リストを与える。2<sup>3</sup>のように通常は同じ素数を複数含むが、
     * 多重度は考慮せずにひとつの素数を１回だけ含むリストを与えればよい。
     * prime_factors が正しくないと結果も正しくないであろう。
     *
     * @param[in] poly GF(2)係数多項式
     * @param[in] degree polyに期待する次数
     * @param[in] prime_factors 2<sup>degree</sup>-1 の素因数分解に現れる素数のリスト
     * @return true poly が原始多項式の場合
     */
    bool isPrime(const NTL::GF2X& poly,
                 int degree, const char * prime_factors[]);

    /**
     * 指定された次数の原始多項式がpolyの因数分解に含まれているか判定する。
     *
     * この関数は、SFMT, dSFMT の開発で使用される。
     *
     * @param[in, out] poly GF(2)係数多項式
     * @param[in] degree polyに期待する次数、メルセンヌ指数であること
     * @return true poly が指定する次数の原始多項式を含む場合
     */
    bool hasFactorOfDegree(NTL::GF2X& poly, long degree);
}
#endif // MTTOOLBOX_PERIOD_HPP
