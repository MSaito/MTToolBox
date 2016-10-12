#ifndef MTTOOLBOX_ALGORITHM_RECURSION_AND_TEMPERING_HPP
#define MTTOOLBOX_ALGORITHM_RECURSION_AND_TEMPERING_HPP
/**
 * @file AlgorithmRecursionAndTempering.hpp
 *
 *\japanese
 * @brief 状態遷移パラメータとテンパリングパラメータの探索を一度に行う
 *\endjapanese
 *
 *\english
 * <ol>
 * <li> Search parameters of state transion function of pseudo
 * random number generator so that the characteristic polynomial
 * of the funciton will have max degree and will be
 * primitive.</li>
 *<li> Search tempering parameters to improve dimension of
 * equi-distribution of output of pseudo random number
 * generator.</li>
 *</ol>
 *\endenglish
 *
 * @author Mutsuo Saito (Hiroshima University)
 * @author Makoto Matsumoto (Hiroshima University)
 *
 * Copyright (C) 2013, 2016 Mutsuo Saito, Makoto Matsumoto and
 * Hiroshima University.
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
#include <MTToolBox/TemperingCalculatable.hpp>
#include <MTToolBox/AlgorithmEquidistribution.hpp>
#include <MTToolBox/AlgorithmRecursionSearch.hpp>
#include <MTToolBox/AlgorithmTempering.hpp>
#include <MTToolBox/util.hpp>

namespace MTToolBox {
    /**
     * @class AlgorithmRecursionAndTempering
     *\japanese
     * \warning このクラスでは、パラメータ探索をする疑似乱数生成器の状態空間のサイズは
     * メルセンヌ指数であると想定されている。
     *<ol>
     * <li>状態遷移関数の特性多項式が原始多項式となるような疑似乱数生成
     * 器のパラメータを探索する。
     * <li>均等分布次元がよくなるようなテンパリングパラメータを探索する。
     *</ol>
     * @tparam U 疑似乱数生成器の出力の型, 符号なし型であること
     * @tparam V パラメータ生成器の出力の型
     *\endjapanese
     *
     *\english
     * <ol>
     * <li> Search parameters of state transion function of pseudo
     * random number generator so that the characteristic polynomial
     * of the funciton will have max degree and will be
     * primitive.</li>
     *<li> Search tempering parameters to improve dimension of
     * equi-distribution of output of pseudo random number
     * generator.</li>
     *</ol>
     * \warning The size of internal state of the generator whose
     * parameters are searched is supporsed to be Mersenne Exponent.
     *
     * @tparam U Type of output of pseudo random number
     * generator. Should be unsigned number.
     * @tparam V Type of output of parameter generator.
     * \endenglish
     */
    template<typename U, typename V = U>
    class AlgorithmRecursionAndTempering {
    public:
        /**
         *\japanese
         * コンストラクタ
         * @param bg パラメータをランダムサーチするための疑似乱数生成器を指定する。
         * この疑似乱数生成器に限り、GF(2)線形である必要はない。疑似乱数生成器である
         * 必要すらない。
         * @param primitivity 原始性判定アルゴリズム。 デフォルトはメルセンヌ
         * 素数周期用のもの。
         *\endjapanese
         *
         *\english
         * Constructor
         * @param bg a generator used for generating numbers to make
         * parameters. This generator is not need to be GF(2)-linear
         * pseudo random number generator, for example, TinyMTDC in
         * sample directory uses sequential counter.
         * @param primitivity An algorithm to check primitivity.
         * As a default, MersennePrimitivity is selected.
         *\endenglish
         */
        AlgorithmRecursionAndTempering(AbstractGenerator<V>& bg,
            const AlgorithmPrimitivity& primitivity = MersennePrimitivity) {
            baseGenerator = &bg;
            isPrime = &primitivity;
        }
        /**
         *\japanese
         * 状態遷移パラメータとテンパリングパラメータを探索する。
         *
         * @param lg テンパリングパラメータ計算可能な疑似乱数生成器
         * @param st1 テンパリングパラメータ探索アルゴリズム
         * @param st2 テンパリングパラメータ探索アルゴリズム（LSB）
         * @param verbose 余分な情報を出力するフラグ
         * @param os 出力ストリーム
         * @param no_lsb LSBからのテンパリングをしない
         * @return true 原始多項式を発見しテンパリングパラメータを設定した場合
         *\endjapanese
         *
         *\english
         * Search parameters for state transition function and
         * parameters for tempering.
         *
         *\endenglish
         * @param lg GF(2)-linear pseudo random number generator
         * whose parameters are to be searched.
         * @param st1 Algorithm for searching tempering parameters.
         * @param st2 Algorithm for searching tempering parameters from LSB.
         * @param verbose if true redundant messages will be outputed.
         * @param os output stream for redundant messages.
         * @param no_lsb if true, \b st2 will not be used.
         * @return false if no tempering parameters which gives proper
         * state transition function are found.
         */
        bool search(TemperingCalculatable<U, V>& lg,
                    AlgorithmTempering<U, V>& st1,
                    AlgorithmTempering<U, V>& st2,
                    bool verbose = false,
                    std::ostream& os = std::cout,
                    bool no_lsb = false) {
            using namespace NTL;
            using namespace std;

            out = &os;
            int veq[bit_size<U>()];
            AlgorithmRecursionSearch<U, V> search(lg, *baseGenerator, *isPrime);
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
                AlgorithmEquidistribution<U, V> sb(lg, bit_size<U>());
                int delta = sb.get_all_equidist(veq);
                print_kv(veq, mexp, bit_size<U>());
                *out << "delta = " << dec << delta << endl;
            }
            if (! no_lsb) {
                st2(lg, verbose);
                if (verbose) {
                    if (st2.isLSBTempering()) {
                        lg.setReverseOutput();
                    }
                    AlgorithmEquidistribution<U, V> sc(lg, bit_size<U>());
                    delta = sc.get_all_equidist(veq);
                    lg.resetReverseOutput();
                    time_t t = time(NULL);
                    *out << "lsb tempering parameters are found at "
                         << ctime(&t) << endl;
                    print_kv(veq, mexp, bit_size<U>());
                    *out << "lsb delta = " << dec << delta << endl;
                }
            }
            st1(lg, verbose);
            AlgorithmEquidistribution<U, V> sc(lg, bit_size<U>());
            delta = sc.get_all_equidist(veq);
            if (verbose) {
                time_t t = time(NULL);
                *out << "tempering parameters are found at " << ctime(&t)
                     << endl;
                *out << lg.getParamString() << endl;
                print_kv(veq, mexp, bit_size<U>());
                *out << "delta = " << dec << delta << endl;
            }
            return true;
        }

        /**
         *\japanese
         *
         * MSBからの均等分布次元のみを向上させたい場合の探索を行う。
         * 状態遷移関数のパラメータは探索する。
         *
         * @param lg テンパリングパラメータ計算可能な疑似乱数生成器
         * @param st テンパリングパラメータ探索アルゴリズム
         * @param verbose 余分な情報を出力するフラグ
         * @param os 出力ストリーム
         * @return 特性多項式が原始多項式となるような状態遷移パラメータを発見した
         *\endjapanese
         *
         *\english
         *
         * Simple wrapper for users who want to search tempering
         * parameter only from MSB.
         *
         * Parameters for state transition function are searched.
         *
         * @param lg GF(2)-linear pseudo random number generator
         * whose parameters are to be searched.
         * @param st Algorithm for searching tempering parameters.
         * @param verbose if true redundant messages will be outputed.
         * @param os output stream for redundant messages.
         * @return false if no tempering parameters which gives proper
         * state transition function are found.
         *\endenglish
         */
        bool search(TemperingCalculatable<U, V>& lg,
                    AlgorithmTempering<U, V>& st,
                    bool verbose = false,
                    std::ostream& os = std::cout) {
            return search(lg, st, st, verbose, os, true);
        }

        /**
         *\japanese
         * 特性多項式のハミングウェイトを返す
         * @return 特性多項式のハミングウェイト
         *\endjapanese
         *
         *\english
         * Returns Hamming weight of characteristic polynomial of state
         * transition function.
         * @return Hamming weight of characteristic polynomial.
         *\endenglish
         */
        int getWeight() {
            return weight;
        }

        /**
         *\japanese
         * 均等分布次元の理論値との差の総和を返す。
         * @return 均等分布次元の理論値との差の総和
         *\endjapanese
         *
         *\english
         * returns sum of d(v)s, which are difference between k(v)
         * and theoretical upper bound.
         * @return sum of d(v)s.
         *\endenglish
         */
        int getDelta() {
            return delta;
        }

        /**
         *\japanese
         * 状態遷移関数の特性多項式を返す。
         * @return 状態遷移関数の特性多項式
         *\endjapanese
         *
         *\english
         * Returns characteristic polynomial of state transition function.
         * @return characteristic polynomial of state transition function.
         *\endenglish
         */
        const NTL::GF2X& getCharacteristicPolynomial() {
            return poly;
        }
    private:
        int weight;
        int delta;
        NTL::GF2X poly;
        std::ostream * out;
        AbstractGenerator<V> * baseGenerator;
        const AlgorithmPrimitivity *isPrime;
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
