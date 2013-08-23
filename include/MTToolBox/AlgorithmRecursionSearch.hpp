#ifndef MTTOOLBOX_ALGORITHM_RECURSION_SEARCH_HPP
#define MTTOOLBOX_ALGORITHM_RECURSION_SEARCH_HPP
/**
 * @file recursion_search.hpp
 *
 * @brief 状態遷移関数のパラメータを探索する。
 *
 * @author Mutsuo Saito (Hiroshima University)
 * @author Makoto Matsumoto (Hiroshima University)
 *
 * Copyright (C) 2013 Mutsuo Saito, Makoto Matsumoto,
 * and Hiroshima University.
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
     * @class AlgorithmRecursionSearch
     * @brief 状態遷移関数のパラメータを探索する。
     *
     * GF(2)線形疑似乱数生成器の状態遷移関数が既約な特性多項式をもつよ
     * うなパラメータを探索する。
     *
     * コンストラクタに与えるgenerator はこのクラスのstart()メソッドを
     * 呼び出すと変更される。そしてstart()がtrueを返したなら、変更後の
     * generatorはその状態遷移関数が既約な特性多項式を持つように変更さ
     * れているはずである。従って変更後のgeneratorからパラメータを取得
     * することによって状態遷移関数が既約となるパラメータを取得できる。
     *
     * <ol><li>start()メソッドを呼ぶ。</li> <li>start()メソッドがtrueを
     * 返したら、getMinpoly(), getCount()などのメソッドを呼んで情報を取
     * 得する。またgeneratorのgetParamString() からも情報を取得できる。
     * </li></ol>
     *
     * @tparam U 疑似乱数生成器の出力する値の型
     */
    template<typename U>
    class AlgorithmRecursionSearch {
    public:
        /**
         * コンストラクタ
         *
         * bg はgenerator とは異なるインスタンスである必要がある。
         * 通常の場合、bg にはMersenneTwisterのインスタンスを与えるとよい。
         * TinyMT では、
         *
         * @param generator generator はこのクラスのメソッドによって変更される。
         * @param bg パラメータ探索に使用する疑似乱数生成器
         */
        AlgorithmRecursionSearch(RecursionSearchable<U>& generator,
                                 AbstractGenerator<U>& bg) {
            rand = &generator;
            baseGenerator = &bg;
            count = 0;
        }

        /**
         *
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
        const std::string getParamString() {
            return rand->getParamString();
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
