#ifndef MTTOOLBOX_ALGORITHM_REDUCIBLE_EQUIDISTRIBUTION_HPP
#define MTTOOLBOX_ALGORITHM_REDUCIBLE_EQUIDISTRIBUTION_HPP
/**
 * @file AlgorithmReducibleEquidistribution.hpp
 *
 *\japanese
 * @brief 可約ジェネレータの均等分布次元を計算する。
 *
 *\endjapanese
 *
 *\english
 * @brief Calculate the dimension of equidistribution of reducible
 * generators.
 *\endenglish
 *
 * @author Mutsuo Saito (Manieth Corp.)
 * @author Makoto Matsumoto (Hiroshima University)
 *
 * Copyright (C) 2015 Mutsuo Saito, Makoto Matsumoto, Manieth Corp.
 * and Hiroshima University.
 * All rights reserved.
 *
 * The 3-clause BSD License is applied to this software, see
 * LICENSE.txt
 */

#include <MTToolBox/ReducibleGenerator.hpp>
#include <MTToolBox/AlgorithmEquidistribution.hpp>

namespace MTToolBox {
    using namespace std;
    /**
     * @class AlgorithmReducibleEquidistribution
     *\japanese
     * @brief 可約ジェネレータの最悪の場合の均等分布次元を計算する。
     *
     * @tparam U 疑似乱数生成器の出力する値の型、符号なし型であること。
     *\endjapanese
     *
     *\english
     * Calculate dimension of equi-distribution of reducible generator
     * in worst case.
     * @tparam U Type of output of pseudo random number
     * generator. Should be unsigned number.
     *\endenglish
     */
    template<typename U, typename G, typename V = U>
    class AlgorithmReducibleEquidistribution {
    public:
        /**
         *\japanese
         * コンストラクタ
         *
         * @param[in] rg
         * @param[in] irreducibleFactor 特性多項式のメルセンヌ指数次の既約成分
         *\endjapanese
         *
         *\english
         * Constructor
         *
         * @param[in] generator GF(2)-linear generator whose
         * dimention of equi-distribution will be searched.
         * @param[in] irreducibleFactor an irreducible factor with
         * Mersenne Exponent degree of characteristic polynomial of
         * generator.
         *\endenglish
         */
        AlgorithmReducibleEquidistribution(const G& rg,
                                           const NTL::GF2X irreducibleFactor,
                                           int bit_length,
                                           int mexp) {
            G * rand = new G(rg);
            NTL::GF2X poly(0,1);
            calcCharacteristicPolynomial(rand, poly);
            NTL::GF2X quotient = poly / irreducibleFactor;
            annihilate<U>(rand, quotient);
            ae = new AlgorithmEquidistribution<U, V>(*rand, bit_length, mexp);
        }

        /**
         *\japanese
         * デストラクタ
         *\endjapanese
         *
         *\english
         * Destructor
         *\endenglish
         */
        ~AlgorithmReducibleEquidistribution() {
            delete ae;
        }

        /**
         *@copydoc AlgorithmEquidistribution::get_all_equidist
         */
        int get_all_equidist(int veq[]) {
            return ae->get_all_equidist(veq);
        }

        /**
         *@copydoc AlgorithmEquidistribution::get_equidist
         */
        int get_equidist(int * sum_equidist) {
            return ae->get_equidist(sum_equidist);
        }
    private:
        AlgorithmEquidistribution<U, V> *ae;
    };
}
#endif // MTTOOLBOX_ALGORITHM_REDUCIBLE_EQUIDISTRIBUTION_HPP
