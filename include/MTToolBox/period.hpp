#ifndef MTTOOLBOX_PERIOD_HPP
#define MTTOOLBOX_PERIOD_HPP
/**
 * @file period.hpp
 *
 *\japanese
 * @brief 最小多項式の計算と原始性の判定
 *\endjapanese
 *
 *\english
 * @brief This file provides functions calculating minimal polynomials
 * and checking primitivity.
 *\endenglish
 *
 * @author Mutsuo Saito (Hiroshima University)
 * @author Makoto Matsumoto (Hiroshima University)
 *
 * Copyright (C) 2013 Mutsuo Saito, Makoto Matsumoto and
 * Hiroshima University.
 * All rights reserved.
 *
 * The 3-clause BSD License is applied to this software, see
 * LICENSE.txt
 */

#include <stdint.h>
#include <NTL/GF2X.h>
#include <NTL/vector.h>
#include <MTToolBox/AbstractGenerator.hpp>
#include <MTToolBox/util.hpp>

namespace MTToolBox {
    /**
     *\japanese
     * 最小多項式を求める
     *
     * @tparam U 疑似乱数生成器の出力の型
     * @param[out] poly 最小多項式
     * @param[in] generator GF(2)疑似乱数生成器
     * @param[in] pos 出力の下位から何ビット目を見るかを指定する
     * @param[in] stateSize 状態空間のビットサイズ
     *\endjapanese
     *
     *\english
     * Calculate minimal polynomial of output sequence.
     *
     * @tparam U type of output of pseudo random number generator
     * @param[out] poly minimal polynomial
     * @param[in] generator GF(2)-linear pseudo random number generator
     * @param[in] pos specifies how manieth bit from LSB is checked,
     * zero means LSB.
     * @param[in] stateSize bit size of internal state.
     *\endenglish
     */
    template<typename U> void
    minpoly(NTL::GF2X& poly, AbstractGenerator<U>& generator, int pos = 0,
            int stateSize = 0)
    {
        using namespace std;
        using namespace NTL;

        Vec<GF2> v;
        int size;
        if (stateSize <= 0) {
            size = generator.bitSize();
        } else {
            size = stateSize;
        }
        v.SetLength(2 * size);
        for (int i = 0; i < 2 * size; i++) {
//            v[i] = (generator.generate() >> pos) & 1;
            v[i] = getBitOfPos(generator.generate(), pos);
        }
        MinPolySeq(poly, v, size);
    }

    /**
     *\japanese
     * 2<sup>degree</sup> -1 が素数となるかどうかを返す
     *
     * メルセンヌ素数の指数のリストを元に判定しているので完全ではない。
     * 疑似乱数生成器として使用する範囲であればカバーしている。
     * @param[in] degree 判定するべき数
     * @return true 2<sup>degree</sup> -1 が素数の場合
     *\endjapanese
     *
     *\english
     * Returns if 2<sup>degree</sup>-1 is prime number.
     *
     * This is checked by fixed list of exponents of Mersenne Prime,
     * therefore this check is not complete.
     * @param[in] degree number to be checked
     * @return true if 2<sup>degree</sup> -1 is prime number.
     *\endenglish
     */
    bool isMexp(uint32_t degree);

    /**
     *\japanese
     * 既約判定
     *
     * @param[in] poly GF(2)係数多項式
     * @return true poly が既約の場合
     *\endjapanese
     *
     *\english
     * Check if polynomial is irreducible
     *
     * @param[in] poly polynomial whose coefficients are GF(2)
     * @return true if \b poly is irreducible
     *\endenglish
     */
    bool isIrreducible(const NTL::GF2X& poly);

    /**
     *\japanese
     * 原始性判定
     *
     * この原始性判定は簡易版であり、poly の次数がメルセンヌ指数の場合のみ
     * 正しい結果を返す。状態空間のビットサイズがメルセンヌ指数でない場合は、
     * この関数を使うべきではない。
     *
     * @param[in] poly GF(2)係数多項式
     * @return true poly の次数がメルセンヌ指数で、かつpolyが既約のとき
     *\endjapanese
     *
     *\english
     * Check if polynomial is primitive.
     *
     * This check is lazy check. Only if degree of polynomial is
     * Mersenne Exponent, the result is correct.
     * This check should not be used, when bit size of internal
     * state is not Mersenne Exponent.
     *
     * @param[in] poly polynomial over GF(2)
     * @return true if polynomial is primitive
     *\endenglish
     */
    bool isPrime(const NTL::GF2X& poly);

    /**
     *\japanese
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
     *\endjapanese
     *
     *\english
     * Check if polynomial is primitive.
     *
     * This function returns false if degree of polynomial is not \b degree.
     * Users should give a list of primes which appear in integer factorization
     * of 2<sup><b>degree</b></sup>-1.
     *
     * @param[in] poly polynomial over GF(2)
     * @param[in] degree \b poly expected to have \b degree.
     * @param[in] prime_factors a list of primes which appear in integer
     * factorization of 2<sup><b>degree</b></sup>-1.
     *\endenglish
     */
    bool isPrime(const NTL::GF2X& poly, int degree,
                 const NTL::Vec<NTL::ZZ>& prime_factors);

    /**
     *\japanese
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
     *\endjapanese
     *
     *\english
     * Check if polynomial is primitive.
     *
     * This function returns false if degree of polynomial is not \b degree.
     * Users should give a list of primes which app-er in integer factorization
     * of 2<sup><b>degree</b></sup>-1.
     *
     * @param[in] poly polynomial over GF(2)
     * @param[in] degree \b poly expected to have \b degree.
     * @param[in] prime_factors a list of primes which appear in integer
     * factorization of 2<sup><b>degree</b></sup>-1.
     *\endenglish
     */
    bool isPrime(const NTL::GF2X& poly,
                 int degree, const char * prime_factors[]);

    /**
     *\japanese
     * 指定された次数の原始多項式がpolyの因数分解に含まれているか判定する。
     *
     * この関数は、SFMT, dSFMT の開発で使用される。
     *
     * @param[in, out] poly GF(2)係数多項式
     * @param[in] degree \b poly に \b degree 次の原始多項式が含まれているか
     * @return true poly が指定する次数の原始多項式を含む場合
     *\endjapanese
     *
     *\english
     * Check if primitive polynomial of given degree appears in
     * factorization of poly.
     *
     * This function will be used by SFMT and dSFMT.
     * @param[in, out] poly polynomial over GF(2)
     * @param[in] degree if \b poly has primitive polynomial with \b degree
     * in its factors.
     * @return true if \b poly has primitive polynomial with \b degree.
     *\endenglish
     */
    bool hasFactorOfDegree(NTL::GF2X& poly, long degree);
}
#endif // MTTOOLBOX_PERIOD_HPP
