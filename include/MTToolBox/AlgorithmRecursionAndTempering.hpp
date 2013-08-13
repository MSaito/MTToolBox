#ifndef MTTOOLBOX_ALGORITHM_RECURSION_AND_TEMPERING_HPP
#define MTTOOLBOX_ALGORITHM_RECURSION_AND_TEMPERING_HPP
/**
 * @file AlgorithmRecursionAndTempering.hpp
 *
 * @brief 状態遷移パラメータとテンパリングパラメータの探索を一度に行う
 *
 * @author Mutsuo Saito (Hiroshima University)
 * @author Makoto Matsumoto (Hiroshima University)
 *
 * Copyright (C) 2013 Mutsuo Saito, Makoto Matsumoto,
 * Hiroshima University and The University of Tokyo.
 * All rights reserved.
 *
 * The 3-clause BSD License is applied to this software, see
 * LICENSE.txt
 */

#include <iostream>
#include <ostream>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <cerrno>
#include <unistd.h>
#include <time.h>
#include <tr1/memory>
#include <MTToolBox/TemperingCalculatable.hpp>
#include <MTToolBox/AlgorithmEquidistribution.hpp>
#include <MTToolBox/AlgorithmRecursionSearch.hpp>
#include <MTToolBox/AlgorithmTempering.hpp>
#include <MTToolBox/util.hpp>

/**
 * @namespace MTToolBox
 * name space for random number generator development
 */
namespace MTToolBox {
    /**
     * @class AlgorithmRecursionAndTempering
     * - search parameters of random number generator whose state transition
     * function has an irreducible characteristic polynomial.
     * - search tempering parameters for the generator.
     *
     * @tparam T type of generator's output, uint32_t or uint64_t.
     * @tparam G class of generator, always linear_generator with template.
     * @tparam ST tempering parameter searching strategy class.
     * @tparam STLSB tempering paramete searching strategy class for
     * tempering from LSB.
     * @tparam SG the sequential generator class
     */
    template<typename T>
    class AlgorithmRecursionAndTempering {
    public:
//        AlgorithmRecursionAndTempering() {
//            const TemperingCalculatable<T>& rand) {
//            generator = rand.copy();
//        }
        /**
         * search
         *
         * @param lg linear generator class
         * @param st tempering parameter searching strategy
         * @param stlsb tempering parameter searching strategy for LSB.
         * @param verbose verbose mode, output information about search process.
         */
        bool search(TemperingCalculatable<T>& lg,
                    AlgorithmTempering<T>& st1,
                    AlgorithmTempering<T>& st2,
                    bool verbose = false,
                    std::ostream& os = std::cout,
                    bool no_lsb = false) {
            using namespace NTL;
            using namespace std;
            using namespace std::tr1;

            out = &os;
            int veq[bit_size(T)];
            AlgorithmRecursionSearch<T> search(lg);
            int mexp = lg.bitSize();
            bool found = false;
            for (int i = 0;; i++) {
                if (search.start(1000 * mexp)) {
                    found = true;
                    break;
                }
                if (verbose) {
                    *out << "not found in " << (i + 1) * 10000 << endl;
                }
            }
            if (!found) {
                return false;
            }
            if (verbose) {
                time_t t = time(NULL);
                *out << "irreducible parameter is found at " << ctime(&t);
            }
            if (verbose) {
                *out << "count = " << search.getCount() << endl;
                lg.printParam(*out);
            }
            poly = search.getMinPoly();
            weight = NTL::weight(poly);
            if (verbose) {
                AlgorithmEquidsitribution<T> sb(lg, bit_size(T));
                //print_binary(*out, poly);
                int delta = sb.get_all_equidist(veq);
                print_kv(veq, mexp, bit_size(T));
                *out << "delta = " << dec << delta << endl;
            }
            if (! no_lsb) {
                st2(lg, verbose);
                if (verbose) {
                    if (st2.isLSBTempering()) {
                        lg.setReverseOutput();
                    }
                    AlgorithmEquidsitribution<T> sc(lg, bit_size(T));
                    delta = sc.get_all_equidist(veq);
                    lg.resetReverseOutput();
                    time_t t = time(NULL);
                    *out << "lsb tempering parameters are found at "
                         << ctime(&t) << endl;
                    print_kv(veq, mexp, bit_size(T));
                    *out << "lsb delta = " << dec << delta << endl;
                }
            }
            st1(lg, verbose);
            AlgorithmEquidsitribution<T> sc(lg, bit_size(T));
            delta = sc.get_all_equidist(veq);
            if (verbose) {
                time_t t = time(NULL);
                *out << "tempering parameters are found at " << ctime(&t)
                     << endl;
                lg.printParam(*out);
                print_kv(veq, mexp, bit_size(T));
                *out << "delta = " << dec << delta << endl;
            }
            return true;
        }

        /**
         * MSBからの均等分布次元のみを向上させたい場合の
         *
         */
        bool search(TemperingCalculatable<T>& lg,
                    AlgorithmTempering<T>& st,
                    bool verbose = false,
                    std::ostream& os = std::cout) {
            return search(lg, st, st, verbose, os, true);
        }

#if 0
        /** getter of rand */
        std::tr1::shared_ptr<TemperingCalculatable<T> > * getGenerator() {
            return generator.copy();
        }
#endif
        /** getter of weight */
        int getWeight() {
            return weight;
        }

        /** getter of delta */
        int getDelta() {
            return delta;
        }

        /** getter of the characteristic polynomial */
        const NTL::GF2X& getCharacteristicPolynomial() {
            return poly;
        }

        /** print the detail of v-bit accuracy equidistribution */
        void print_kv(int veq[], int mexp, int size) {
            using namespace std;
            for (int i = 0; i < size; i++) {
                *out << dec << i + 1 << ":" << veq[i]
                     << "(" << mexp / (i + 1) - veq[i] << ")"
                     << endl;
            }
        }
    private:
        /**
         * the humming weight of the charateristic polynomial.
         * e.g the number of terms of characteristic polynomial.
         */
        int weight;
        /**
         * sum of d(v) for all \b v.
         */
        int delta;
        /**
         * searched generator with its parameters.
         */
        //std::tr1::shared_ptr<TemperingCalculatable<T> > generator;
        /**
         * the characteristic polynomial of linear transition of the
         * internal state.
         */
        NTL::GF2X poly;
        std::ostream *out;
    };

}
#endif //MTTOOLBOX_ALGORITHM_RECURSION_AND_TEMPERING_HPP
