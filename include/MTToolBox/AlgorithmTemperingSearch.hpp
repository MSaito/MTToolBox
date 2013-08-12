#ifndef MTTOOLBOX_ALGORITHM_TEMPERING_SEARCH_HPP
#define MTTOOLBOX_ALGORITHM_TEMPERING_SEARCH_HPP
/**
 * @file AlgorithmTemperingSearch.hpp
 *
 * @brief テンパリングパラメータを探索する
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

#include <ostream>
#include <MTToolBox/TemperingSearchable.hpp>
#include <MTToolBox/AlgothmTempering.hpp>


namespace MTToolBox {
    /**
     * @class Search
     * search parameters so that the generator's state transition function
     * has an irreducible characteristic polynomial.
     * 1) call start() function.
     * 2) if start() returns true, then call get_random(), get_minpoly(),
     * or get_count().
     *
     * @tparam T generators class
     */
    template<typename U>
    class AlgorithmTemperingSearch {
    public:
        /**
         * rand は破壊される。
         * @param rand_ random number generator whose parameters are
         * searched.
         * @param seq_generator a sequential number generator which
         * gives sequential number for searching parameters.
         *
         */
        AlgorithmTemperingSearch(TemperingSearchable<U>& rand_) {
            rand = &rand_;
        }

        /**
         * generate random parameters and check if the generator's
         * state transition function has an irreducible characteristic
         * polynomial. If found in \b try_count times, return true, else
         * return false.
         *
         * @param try_count
         */
        bool operator()(AlgorithmTempering& alg) {
            long size = rand->bitSize();
            long degree;
            for (int i = 0; i < try_count; i++) {
                rand->setUpParam();
                rand->seed(1);
                minpoly(poly, *rand);
                count++;
                degree = deg(poly);
                if (degree != size) {
                    continue;
                }
                if (isIrreducible(poly)) {
                    return true;
                }
            }
            return false;
        }

    private:
        RecursionSearchable<U> *rand;
        NTL::GF2X poly;
        long count;
    };
}
#endif // MTTOOLBOX_ALGORITHM_TEMPERING_SEARCH_HPP
