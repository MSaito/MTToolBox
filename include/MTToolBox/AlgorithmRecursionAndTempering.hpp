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

namespace MTToolBox {
    /**
     * @class AlgorithmRecursionAndTempering
     * - 状態遷移関数の特性多項式が原始多項式となるような疑似乱数生成器のパラメータを探索する。
     * - 均等分布次元がよくなるようなテンパリングパラメータを探索する。
     *
     * @tparam T 疑似乱数生成器の出力の型
     */
    template<typename T>
    class AlgorithmRecursionAndTempering {
    public:
        AlgorithmRecursionAndTempering(AbstractGenerator<T>& bg) {
            baseGenerator = &bg;
        }
        /**
         * 状態遷移パラメータとテンパリングパラメータを探索する。
         *
         * @param lg テンパリングパラメータ計算可能な疑似乱数生成器
         * @param st1 テンパリングパラメータ探索アルゴリズム
         * @param st2 テンパリングパラメータ探索アルゴリズム（LSB）
         * @param verbose 余分な情報を出力するフラグ
         * @param os 出力ストリーム
         * @param no_lsb LSBからのテンパリングをしない
         * @return 原始多項式を発見した
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
            int veq[bit_size<T>()];
            AlgorithmRecursionSearch<T> search(lg, *baseGenerator);
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
                *out << lg.getParamString() << endl;
            }
            poly = search.getMinPoly();
            weight = NTL::weight(poly);
            if (verbose) {
                AlgorithmEquidsitribution<T> sb(lg, bit_size<T>());
                //print_binary(*out, poly);
                int delta = sb.get_all_equidist(veq);
                print_kv(veq, mexp, bit_size<T>());
                *out << "delta = " << dec << delta << endl;
            }
            if (! no_lsb) {
                st2(lg, verbose);
                if (verbose) {
                    if (st2.isLSBTempering()) {
                        lg.setReverseOutput();
                    }
                    AlgorithmEquidsitribution<T> sc(lg, bit_size<T>());
                    delta = sc.get_all_equidist(veq);
                    lg.resetReverseOutput();
                    time_t t = time(NULL);
                    *out << "lsb tempering parameters are found at "
                         << ctime(&t) << endl;
                    print_kv(veq, mexp, bit_size<T>());
                    *out << "lsb delta = " << dec << delta << endl;
                }
            }
            st1(lg, verbose);
            AlgorithmEquidsitribution<T> sc(lg, bit_size<T>());
            delta = sc.get_all_equidist(veq);
            if (verbose) {
                time_t t = time(NULL);
                *out << "tempering parameters are found at " << ctime(&t)
                     << endl;
                *out << lg.getParamString() << endl;
                print_kv(veq, mexp, bit_size<T>());
                *out << "delta = " << dec << delta << endl;
            }
            return true;
        }

        /**
         * MSBからの均等分布次元のみを向上させたい場合の探索を行う。
         *
         * @param lg テンパリングパラメータ計算可能な疑似乱数生成器
         * @param st テンパリングパラメータ探索アルゴリズム
         * @param verbose 余分な情報を出力するフラグ
         * @param os 出力ストリーム
         * @return 特性多項式が原始多項式となるような状態遷移パラメータを発見した
         */
        bool search(TemperingCalculatable<T>& lg,
                    AlgorithmTempering<T>& st,
                    bool verbose = false,
                    std::ostream& os = std::cout) {
            return search(lg, st, st, verbose, os, true);
        }

        /** 特性多項式のハミングウェイトを返す */
        int getWeight() {
            return weight;
        }

        /** 均等分布次元の理論値との差の総和を返す。*/
        int getDelta() {
            return delta;
        }

        /** 状態遷移関数の特性多項式を返す。*/
        const NTL::GF2X& getCharacteristicPolynomial() {
            return poly;
        }
    private:
        int weight;
        int delta;
        NTL::GF2X poly;
        std::ostream * out;
        AbstractGenerator<T> * baseGenerator;

        void print_kv(int veq[], int mexp, int size) {
            using namespace std;
            for (int i = 0; i < size; i++) {
                *out << dec << i + 1 << ":" << veq[i]
                     << "(" << mexp / (i + 1) - veq[i] << ")"
                     << endl;
            }
        }
    };

}
#endif //MTTOOLBOX_ALGORITHM_RECURSION_AND_TEMPERING_HPP
