#ifndef MTTOOLBOX_ALGORITHM_RECURSION_SEARCH_HPP
#define MTTOOLBOX_ALGORITHM_RECURSION_SEARCH_HPP
/**
 * @file recursion_search.hpp
 *
 * @brief search parameters so that the random number generator's state
 * transition function has an irreducible characteristic polynomial.
 *
 * @author Mutsuo Saito (Hiroshima University)
 * @author Makoto Matsumoto (The University of Tokyo)
 *
 * Copyright (C) 2011 Mutsuo Saito, Makoto Matsumoto,
 * Hiroshima University and The University of Tokyo.
 * All rights reserved.
 *
 * The 3-clause BSD License is applied to this software, see
 * LICENSE.txt
 */

#include <ostream>
#include <NTL/GF2X.h>
#include <NTL/GF2XFactoring.h>
#include <MTToolBox/RecursionSearchable.hpp>
#include <MTToolBox/period.hpp>

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
    class AlgorithmRecursionSearch {
    public:
        /**
         * rand は破壊される。
         * @param rand_ random number generator whose parameters are
         * searched.
         * @param seq_generator a sequential number generator which
         * gives sequential number for searching parameters.
         *
         */
        AlgorithmRecursionSearch(RecursionSearchable<U>& rand_,
                                 AbstractGenerator<U>& bg) {
            rand = &rand_;
            baseGenerator = &bg;
            count = 0;
        }

        /**
         * generate random parameters and check if the generator's
         * state transition function has an irreducible characteristic
         * polynomial. If found in \b try_count times, return true, else
         * return false.
         *
         * @param try_count
         */
        bool start(int try_count) {
            long size = rand->bitSize();
            long degree;
            for (int i = 0; i < try_count; i++) {
                rand->setUpParam(*baseGenerator);
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

        /**
         * call this function after \b start() has returned true.
         * @return random number generator class with parameters.
         */
        void printParam(std::ostream& out) {
            rand->printParam(out);
        }

        /**
         * call this function after \b start() has returned true.
         * In this program, if minimal polynomial is irreducible,
         * then the polynomial is characteristic polynomial of
         * generator's state transition function.
         *
         * @return minimal polynomial of generated sequence.
         */
        const NTL::GF2X& getMinPoly() const {
            return poly;
        }

        /**
         * @return tried count after this class has created.
         */
        long getCount() const {
            return count;
        }

    private:
        RecursionSearchable<U> *rand;
        AbstractGenerator<U> *baseGenerator;
        NTL::GF2X poly;
        long count;
    };
}
#endif // MTTOOLBOX_ALGORITHM_RECURSION_SEARCH_HPP
