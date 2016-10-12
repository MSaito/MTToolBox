#ifndef MTTOOLBOX_ALGORITHM_REDUCIBLE_RECURSION_SEARCH_HPP
#define MTTOOLBOX_ALGORITHM_REDUCIBLE_RECURSION_SEARCH_HPP
/**
 * @file AlgorithmReducibleRecursionSearch.hpp
 *
 *\japanese
 * @brief 可約ジェネレータの状態遷移関数のパラメータを探索する。
 *
 * F_2疑似乱数生成器の状態遷移関数の特性多項式が、「指定されたメルセンヌ
 * 指数次でかつ規約な因子」を持つようなパラメータを探索する。
 *
 *\endjapanese
 *
 *\english
 * @brief Search parameters of state transition function.  Search
 * parameters of state transition function of reducible pseudo random
 * number generator so that the characteristic polynomial of the
 * function will have a factor polynomial which has specified mersenne
 * expornent degree and is irreducible.
 *\endenglish
 *
 * @author Mutsuo Saito (Manieth Corp.)
 * @author Makoto Matsumoto (Hiroshima University)
 *
 * Copyright (C) 2014, 2016 Mutsuo Saito, Makoto Matsumoto, Manieth Corp.
 * and Hiroshima University.
 * All rights reserved.
 *
 * The 3-clause BSD License is applied to this software, see
 * LICENSE.txt
 */

#if defined(DEBUG)
#include <iostream>
#endif
#include <NTL/GF2X.h>
#include <NTL/GF2XFactoring.h>
#include <MTToolBox/util.hpp>
#include <MTToolBox/ReducibleGenerator.hpp>
#include <MTToolBox/AlgorithmPrimitivity.hpp>
#include <MTToolBox/period.hpp>

namespace MTToolBox {
    using namespace std;

    template<typename U, typename V = U>
    void calcCharacteristicPolynomial(RecursionSearchable<U, V> *rand,
                                      NTL::GF2X& poly);
    /**
     * @class AlgorithmReducibleRecursionSearch
     *\japanese
     * @brief 可約ジェネレータの状態遷移関数のパラメータを探索する。
     *
     * GF(2)線形疑似乱数生成器の状態遷移関数の特性多項式が既約なメルセ
     * ンヌ指数次の因子をもつをもつようなパラメータを探索する。
     *
     * コンストラクタに与えるgenerator はこのクラスのstart()メソッドを
     * 呼び出すと変更される。そしてstart()がtrueを返したなら、変更後の
     * generatorはその状態遷移関数が既約な特性多項式を持つように変更さ
     * れているはずである。従って変更後のgeneratorからパラメータを取得
     * することによって状態遷移関数が既約となるパラメータを取得できる。
     *
     * <ol><li>start()メソッドを呼ぶ。</li> <li>start()メソッドがtrueを
     * 返したら、getMinpoly(), getIrreducibleFactor(),
     * getCount()などのメソッドを呼んで情報を取
     * 得する。またgeneratorのgetParamString() からも情報を取得できる。
     * </li></ol>
     *
     * @tparam U 疑似乱数生成器の出力する値の型、符号なし型であること。
     * @tparam V パラメータ生成用疑似乱数生成器の出力する値の型
     *\endjapanese
     *
     *\english
     * Search parameters of state transition function of pseudo
     * random number generator so that the characteristic polynomial
     * of the function will have max degree and will be
     * primitive.
     * @tparam U Type of output of pseudo random number
     * generator. Should be unsigned number.
     * @tparam V Type of output of base generator. Should be unsigned number.
     *\endenglish
     */
    template<typename U, typename V = U>
    class AlgorithmReducibleRecursionSearch {
    public:
        /**
         *\japanese
         * コンストラクタ
         *
         * このコンストラクタでは、状態空間の大きさがメルセンヌ指数より
         * 少し大きい２のべき乗の場合を想定している。
         *
         * bg はgenerator とは異なるインスタンスである必要がある。
         * 通常の場合、bg にはMersenneTwisterのインスタンスを与えるとよい。
         * TinyMT では、規則的に減少する数列を使用している。
         *
         * @param generator generator はこのクラスのメソッドによって変更される。
         * @param bg パラメータ探索に使用する疑似乱数生成器
         *\endjapanese
         *
         *\english
         * Constructor
         * Size of internal state is supposed to be power of 2 and
         * larger than Mersenne Exponent.
         * @param[in,out] generator GF(2)-linear generator whose parameters will
         * be searched. The generator will be changed.
         * @param[in,out] bg a generator used for generating numbers to make
         * parameters. This generator is not need to be GF(2)-linear
         * pseudo random number generator.
         *\endenglish
         */
        AlgorithmReducibleRecursionSearch(ReducibleGenerator<U, V>& generator,
                                 AbstractGenerator<V>& bg) {
            rand = &generator;
            baseGenerator = &bg;
            count = 0;
        }

        /**
         * @copydoc AlgorithmRecursionSearch::start
         */
        bool start(int try_count) {
            long degree;
            long mexp = rand->getMexp();
            for (int i = 0; i < try_count; i++) {
                rand->setUpParam(*baseGenerator);
                rand->seed(getOne<U>());
#if defined(DEBUG)
                cout << "rand param:";
                cout << rand->getParamString() << endl;
#endif
                minpoly(poly, *rand);
                irreducible = poly;
                count++;
                if (deg(irreducible) < mexp) {
#if defined(DEBUG)
                    cout << "irrepoly deg = " << deg(irreducible)
                         << " skip" << endl;
#endif
                    continue;
                }
                bool hasFactor = hasFactorOfDegree(irreducible, mexp);
#if defined(DEBUG)
                cout << "has factor = " << hasFactor << endl;
#endif
                if (!hasFactor) {
#if defined(DEBUG)
                    cout << "not has factor of degree irrepoly deg = "
                         << deg(irreducible) << " skip" << endl;
#endif
                    continue;
                }
                degree = deg(irreducible);
                if (degree != mexp) {
#if defined(DEBUG)
                    cout << "degree:" << degree << "degree != mexp skip"
                         << endl;
#endif
                    continue;
                }
                return true;
            }
            return false;
        }

        /**
         * @copydoc AlgorithmRecursionSearch::getParamString
         */
        const std::string getParamString() {
            return rand->getParamString();
        }

        /**
         * @copydoc AlgorithmRecursionSearch::getMinPoly
         *\japanese
         * このメソッドが返すのは特性多項式ではない可能性がある。
         *\endjapanese
         *
         *\english
         * This method may return a polynomial which may not be a
         * characteristic polynomial.
         *\endenglish
         */
        const NTL::GF2X& getCharacteristicPolynomial() {
            return poly;
        }

        /**
         *\japanese
         * 特性多項式のメルセンヌ指数次の大きな既約因子を返す
         * このメソッドは start() が true を返した場合にのみ呼び出すべきである。
         *
         * @return 最大既約因子
         *\endjapanese
         *
         *\english
         * Returns a maxmum irreducible factor of minimal polynomial
         * of output of the pseudo random number generator.
         * Call this method only after start() returns true.
         * @return a maxmun irreducible factor of minimal polynomial
         *\endenglish
         */
        const NTL::GF2X& getIrreducibleFactor() const {
            return irreducible;
        }

        /**
         *\japanese
         * このインスタンスが作られてからstart() が終了するまでに試行した回数を返す。
         * @return このインスタンスが作られてからstart() が終了するまでに試行した回数
         *\endjapanese
         *
         *\english
         * Returns tried count from the instance was created.
         * @return tried count from the instance was created.
         *\endenglish
         */
        long getCount() const {
            return count;
        }

    private:
        ReducibleGenerator<U, V> *rand;
        AbstractGenerator<V> *baseGenerator;
        NTL::GF2X poly;
        NTL::GF2X irreducible;
        long count;

    };

    /**
     *\japanese
     * Reducible Generator の特性多項式の計算
     * 実のところ、特性多項式ではなく最小多項式のLCMを計算しているに過ぎない。
     * 次数が一致すれば特性多項式。特性多項式でなくても、MTToolBoxで使用する
     * 範囲内では特に問題はない。
     * @tparam U 疑似乱数生成器の返す値の型
     * @tparam V パラメータ生成器の返す値の型
     * @param[in,out] rand 疑似乱数生成器
     * @param[in,out] poly 特性多項式
     *\endjapanese
     *
     *\english
     * Calculate the characteristic polynomial of reducible generator.
     * @tparam U type of return value of random number generator.
     * @tparam V type of return value of parameter generator.
     * @param[in,out] rand pseudo random number generator.
     * @param[in,out] poly calculated characteristic polynomial.
     *\endenglish
     */
    template<typename U, typename V = U>
    void calcCharacteristicPolynomial(ReducibleGenerator<U, V> *rand,
                                      NTL::GF2X& poly)
    {
#if defined(DEBUG)
        cout << "calcCharacteristicPolynomial start" << endl;
#endif
        int bitsize = bit_size<U>();
        NTL::GF2X minpol;
        NTL::GF2X lcmpoly = poly;
        int size = rand->bitSize();
        for (int i = 0; i < bitsize; i++) {
            if (deg(lcmpoly) == size) {
                poly = lcmpoly;
#if defined(DEBUG)
                cout << "calcCharacteristicPolynomial end" << endl;
#endif
                return;
            }
            minpoly(minpol, *rand, i);
            LCM(lcmpoly, lcmpoly, minpol);
        }
        poly = lcmpoly;
#if defined(DEBUG)
        cout << "calcCharacteristicPolynomial end" << endl;
#endif
    }
}
#endif // MTTOOLBOX_ALGORITHM_REDUCIBLE_RECURSION_SEARCH_HPP
