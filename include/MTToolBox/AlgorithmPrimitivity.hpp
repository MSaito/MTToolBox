#ifndef MTTOOLBOX_ALGORITHM_PRIMITIVITY_HPP
#define MTTOOLBOX_ALGORITHM_PRIMITIVITY_HPP
/**
 * @file AlgorithmPrimitivity.hpp
 *
 * @brief 原始多項式判定
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
     * @brief 原始多項式かどうか判定するアルゴリズムを提供するクラス
     *
     */
    class AlgorithmPrimitivity {
    public:
        /**
         * メルセンヌ指数を次数とする原始多項式かどうか判定する場合のコンストラクタ
         *
         * operator() によって多項式は degree 次数でかつ既約かテストされる。
         */
        AlgorithmPrimitivity() {
            primes = new NTL::Vec<NTL::ZZ>;
            primes->SetLength(0);
            mersenne = true;
        }

        /**
         * 一般のGF(2)係数の原始多項式かどうか判定する場合のコンストラクタ
         * @param[in] prime_factors 2<sup>degree</sup> -1の素因数分解に現れる
         * 素数の文字列表現のリスト
         */
        AlgorithmPrimitivity(const char * prime_factors[]);

        /**
         * デストラクタ
         */
        ~AlgorithmPrimitivity() {
            delete primes;
        }

        /**
         * poly がコンストラクタで指定した次数の原始多項式かどうか判定する
         *
         * @param[in] max_degree 状態空間の大きさから定まる最大次数
         * @param[in] poly GF(2)係数多項式
         * @return true 最大次数の原始多項式の場合
         */
        bool operator()(int max_degree, const NTL::GF2X& poly) const;
    private:
        bool mersenne;
        NTL::Vec<NTL::ZZ> * primes;
    };

    extern const AlgorithmPrimitivity MersennePrimitivity;
}
#endif // MTTOOLBOX_ALGORITHM_PRIMITIVITY_HPP
