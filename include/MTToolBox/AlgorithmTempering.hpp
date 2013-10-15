#ifndef MTTOOLBOX_ALGORITHM_TEMPERING_HPP
#define MTTOOLBOX_ALGORITHM_TEMPERING_HPP
/**
 * @file AlgorithmTempering.hpp
 *
 *\japanese
 * @brief テンパリングパラメータ探索アルゴリズムのための抽象クラス
 *\endjapanese
 *
 *\english
 * @brief Abstract class for searching tempering parameters.
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
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <unistd.h>
#include <tr1/memory>
#include <MTToolBox/TemperingCalculatable.hpp>

namespace MTToolBox {
    /**
     * @class AlgorithmTempering
     *
     *\japanese
     * @brief 疑似乱数生成器の高次元均等分布性を改善するために、テンパ
     * リングパラメータを探索するアルゴリズム
     *
     * \warning テンパリングパラメータの探索をしても十分良い高次元均等
     * 分布が得られない場合は、状態遷移関数の変更を考慮した方がよい。状
     * 態遷移関数で十分ビットミックスされていない場合、単純なテンパリン
     * グで均等分布次元を最大化することはできないだろう。
     *
     * @tparam U 疑似乱数生成器の出力の型, 例えば uint32_t など。
     *\endjapanese
     *
     *\english
     * @brief Algorithm that search tempering parameters to improve
     * dimension of equi-distribution of output of pseudo random number
     * generator.
     * \warning If you could not get high dimension of
     * equi-distribution after tempering using parameters got by this
     * algorithm, you should consider changing the design of random
     * number generation algorithm.
     *
     * @tparam U type of output of pseudo random number generator,
     * for example, uint32_t. Only unsigned numbers are allowed.
     *\endenglish
     */
    template<typename U>
    class AlgorithmTempering {
    public:

        /**
         *\japanese
         * テンパリングパラメータを探索する
         * @param[in, out] rand 疑似乱数生成器
         * @param[in] verbose 余分な情報を表示する
         * @returns 0
         *\endjapanese
         *
         *\english
         * Search tempering parameters. Searched parameters are
         * set to \b rand.
         * This process may consume large CPU time.
         * @param rand pseudo random number generator
         * @param verbose if true output redundant messages.
         * @return always zero.
         *\endenglish
         */
        virtual int operator()(TemperingCalculatable<U>& rand,
                               bool verbose = false) = 0;

        /**
         *\japanese
         * LSB からのテンパリングをするのか
         * @return true LSBからのテンパリング
         *\endjapanese
         *
         *\english
         * Shows if searching tempering parameters is from LSBs.
         * @return true if searching is from LSBs.
         *\endenglish
         */
        virtual bool isLSBTempering() const {
            return false;
        }
    };
}

#endif // MTTOOLBOX_ALGORITHM_TEMPERING_HPP
