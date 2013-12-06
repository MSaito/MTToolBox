#ifndef MTTOOLBOX_ALGORITHM_PRIMITIVITY_HPP
#define MTTOOLBOX_ALGORITHM_PRIMITIVITY_HPP
/**
 * @file AlgorithmPrimitivity.hpp
 *
 *\japanese
 * @brief 原始多項式判定
 *\endjapanese
 *
 *\english
 * @brief Primitivity test
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
 */

#include <NTL/GF2X.h>
#include <NTL/ZZ.h>
#include <NTL/vector.h>
#include <MTToolBox/period.hpp>

namespace MTToolBox {
    /**
     * @class AlgorithmPrimitivity
     *
     *\japanese
     * @brief 原始多項式かどうか判定するアルゴリズムを提供するクラス
     *\endjapanese
     *
     *\english
     * Algorithm which check if given polynomial is primitive.
     *\endenglish
     */
    class AlgorithmPrimitivity {
    public:
        /**
         *\japanese
         * メルセンヌ指数を次数とする原始多項式かどうか判定する場合のコンストラクタ
         *
         * operator() によって多項式は degree 次数でかつ既約かテストされる。
         *\endjapanese
         *
         *\english
         * Constructor for polynomials whose degrees are Mersenne Exponent.
         *\endenglish
         */
        AlgorithmPrimitivity() {
            primes = new NTL::Vec<NTL::ZZ>;
            primes->SetLength(0);
            mersenne = true;
        }

        /**
         *\japanese
         * 一般のGF(2)係数の原始多項式かどうか判定する場合のコンストラクタ
         * @param[in] prime_factors 2<sup>degree</sup> -1の素因数分解に現れる
         * 素数の文字列表現のリスト
         *\endjapanese
         *
         *\english
         * Constructor for general polynomial with GF(2) coefficients.
         * @param[in] prime_factors List of character strings of primes
         * which appear in integer factorization of 2<sup>degree</sup>-1.
         *\endenglish
         */
        AlgorithmPrimitivity(const char * prime_factors[]);

        /**
         *\japanese
         * デストラクタ
         *\endjapanese
         *\english
         * Destructor
         *\endenglish
         */
        ~AlgorithmPrimitivity() {
            delete primes;
        }

        /**
         *\japanese
         * poly が指定した次数の原始多項式かどうか判定する
         *
         * @param[in] max_degree 状態空間の大きさから定まる最大次数
         * @param[in] poly GF(2)係数多項式
         * @return true 最大次数の原始多項式の場合
         *\endjapanese
         *
         *\english
         * Check if given polynomial is a primitive polynomial of given \b
         * max_degree.
         *
         * @param[in] max_degree 状態空間の大きさから定まる最大次数
         * @param[in] poly GF(2)係数多項式
         * @return true 最大次数の原始多項式の場合
         *\endenglish
         */
        bool operator()(int max_degree, const NTL::GF2X& poly) const;
    private:
        bool mersenne;
        NTL::Vec<NTL::ZZ> * primes;
    };

    /**
     *\japanese
     * 状態空間のビットサイズがメルセンヌ指数の疑似乱数生成器の
     * 最小多項式の原始性を判定するアルゴリズム
     *\endjapanese
     *
     *\english
     * An algorithm which checks if given polynomial is a
     * primitive polynomial of given degree for pseudo
     * random number generator whose internal state size
     * is Mersenne exponent.
     *\endenglish
     */
    extern const AlgorithmPrimitivity MersennePrimitivity;

    /**
     *\japanese
     * 2<sup>128</sup>-1 の素因数分解に現れる素数のリスト
     *\endjapanese
     *
     *\english
     * List of prime numbers appear in the factorization
     * of 2<sup>128</sup>-1.
     *\endenglish
     */
    extern const char * prime_factors2_128_1[];

    /**
     *\japanese
     * 2<sup>160</sup>-1 の素因数分解に現れる素数のリスト
     *\endjapanese
     *
     *\english
     * List of prime numbers appear in the factorization
     * of 2<sup>160</sup>-1.
     *\endenglish
     */
    extern const char * prime_factors2_160_1[];

    /**
     *\japanese
     * 2<sup>192</sup>-1 の素因数分解に現れる素数のリスト
     *\endjapanese
     *
     *\english
     * List of prime numbers appear in the factorization
     * of 2<sup>192</sup>-1.
     *\endenglish
     */
    extern const char * prime_factors2_192_1[];

    /**
     *\japanese
     * 2<sup>224</sup>-1 の素因数分解に現れる素数のリスト
     *\endjapanese
     *
     *\english
     * List of prime numbers appear in the factorization
     * of 2<sup>224</sup>-1.
     *\endenglish
     */
    extern const char * prime_factors2_224_1[];

    /**
     *\japanese
     * 2<sup>256</sup>-1 の素因数分解に現れる素数のリスト
     *\endjapanese
     *
     *\english
     * List of prime numbers appear in the factorization
     * of 2<sup>256</sup>-1.
     *\endenglish
     */
    extern const char * prime_factors2_256_1[];

    /**
     *\japanese
     * 2<sup>288</sup>-1 の素因数分解に現れる素数のリスト
     *\endjapanese
     *
     *\english
     * List of prime numbers appear in the factorization
     * of 2<sup>288</sup>-1.
     *\endenglish
     */
    extern const char * prime_factors2_288_1[];

    /**
     *\japanese
     * 2<sup>320</sup>-1 の素因数分解に現れる素数のリスト
     *\endjapanese
     *
     *\english
     * List of prime numbers appear in the factorization
     * of 2<sup>320</sup>-1.
     *\endenglish
     */
    extern const char * prime_factors2_320_1[];

    /**
     *\japanese
     * 2<sup>352</sup>-1 の素因数分解に現れる素数のリスト
     *\endjapanese
     *
     *\english
     * List of prime numbers appear in the factorization
     * of 2<sup>352</sup>-1.
     *\endenglish
     */
    extern const char * prime_factors2_352_1[];

    /**
     *\japanese
     * 2<sup>384</sup>-1 の素因数分解に現れる素数のリスト
     *\endjapanese
     *
     *\english
     * List of prime numbers appear in the factorization
     * of 2<sup>384</sup>-1.
     *\endenglish
     */
    extern const char * prime_factors2_384_1[];

    /**
     *\japanese
     * 2<sup>416</sup>-1 の素因数分解に現れる素数のリスト
     *\endjapanese
     *
     *\english
     * List of prime numbers appear in the factorization
     * of 2<sup>416</sup>-1.
     *\endenglish
     */
    extern const char * prime_factors2_416_1[];

    /**
     *\japanese
     * 2<sup>448</sup>-1 の素因数分解に現れる素数のリスト
     *\endjapanese
     *
     *\english
     * List of prime numbers appear in the factorization
     * of 2<sup>448</sup>-1.
     *\endenglish
     */
    extern const char * prime_factors2_448_1[];

    /**
     *\japanese
     * 2<sup>480</sup>-1 の素因数分解に現れる素数のリスト
     *\endjapanese
     *
     *\english
     * List of prime numbers appear in the factorization
     * of 2<sup>480</sup>-1.
     *\endenglish
     */
    extern const char * prime_factors2_480_1[];

    /**
     *\japanese
     * 2<sup>512</sup>-1 の素因数分解に現れる素数のリスト
     *\endjapanese
     *
     *\english
     * List of prime numbers appear in the factorization
     * of 2<sup>512</sup>-1.
     *\endenglish
     */
    extern const char * prime_factors2_512_1[];

    /**
     *\japanese
     * 2<sup>544</sup>-1 の素因数分解に現れる素数のリスト
     *\endjapanese
     *
     *\english
     * List of prime numbers appear in the factorization
     * of 2<sup>544</sup>-1.
     *\endenglish
     */
    extern const char * prime_factors2_544_1[];

}
#endif // MTTOOLBOX_ALGORITHM_PRIMITIVITY_HPP
