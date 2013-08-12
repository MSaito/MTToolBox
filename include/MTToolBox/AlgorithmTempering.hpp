#ifndef MTTOOLBOX_ALGORITHM_TEMPERING_HPP
#define MTTOOLBOX_ALGORITHM_TEMPERING_HPP
/**
 * @file AlgorithmTempering.hpp
 *
 * @brief search tempering parameters so that the random number
 * generator has a good equidistribution property.
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
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <unistd.h>
#include <tr1/memory>
#include <MTToolBox/abstract_temper_searcher.hpp>
#include <MTToolBox/calc_equidist.hpp>

namespace MTToolBox {
    /**
     * @class AlgorithmTempering
     *
     * @brief 疑似乱数生成器の高次元均等分布性を改善するために、テンパ
     * リングパラメータを探索するアルゴリズム
     *
     * @caution テンパリングパラメータの探索をしても十分良い高次元均等
     * 分布が得られない場合は、状態遷移関数の変更を考慮した方がよい。状
     * 態遷移関数で十分ビットミックスされていない場合、単純なテンパリン
     * グで均等分布次元を最大化することはできないだろう。
     *
     * @tparam U 疑似乱数生成器の出力の型, 例えば uint32_t など。
     */
    template<typename U>
    class AlgorithmTempering {
    public:

        /**
         * search tempering parameters.
         * @returns 0
         */
        int operator()(TemperingSearchable<U>& rand,
                       bool verbose = false);
    };
}

#endif // MTTOOLBOX_ALGORITHM_TEMPERING_HPP
