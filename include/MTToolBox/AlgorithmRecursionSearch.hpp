#ifndef MTTOOLBOX_ALGORITHM_RECURSION_SEARCH_HPP
#define MTTOOLBOX_ALGORITHM_RECURSION_SEARCH_HPP
/**
 * @file AlgorithmRecursionSearch.hpp
 *
 *\japanese
 * @brief 状態遷移関数のパラメータを探索する。
 *\endjapanese
 *
 *\english
 * @brief Search parameters of state transition function.
 * Search parameters of state transition function of pseudo random
 * number generator so that the characteristic polynomial of the
 * function will have max degree and will be primitive.
 *\endenglish
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

#if defined(DEBUG)
#include <iostream>
#endif
#include <NTL/GF2X.h>
#include <NTL/GF2XFactoring.h>
#include <MTToolBox/RecursionSearchable.hpp>
#include <MTToolBox/AlgorithmPrimitivity.hpp>
#include <MTToolBox/period.hpp>
#include <MTToolBox/util.hpp>

namespace MTToolBox {
    using namespace std;
    /**
     * @class AlgorithmRecursionSearch
     *\japanese
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
     * @tparam U 疑似乱数生成器の出力する値の型、符号なし型であること。
     *\endjapanese
     *
     *\english
     * Search parameters of state transition function of pseudo
     * random number generator so that the characteristic polynomial
     * of the function will have max degree and will be
     * primitive.
     * @tparam U Type of output of pseudo random number
     * generator. Should be unsigned number.
     *\endenglish
     */
    template<typename U, typename V = U>
    class AlgorithmRecursionSearch {
    public:
        /**
         *\japanese
         * コンストラクタ
         *
         * このコンストラクタでは、状態空間の大きさがメルセンヌ指数の場合を
         * 想定している。
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
         * Size of internal state is supposed to be Mersenne Exponent.
         * @param[in,out] generator GF(2)-linear generator whose parameters will
         * be searched. The generator will be changed.
         * @param[in,out] bg a generator used for generating numbers to make
         * parameters. This generator is not need to be GF(2)-linear
         * pseudo random number generator, for example, TinyMTDC in
         * sample directory uses sequential counter.
         *\endenglish
         */
        AlgorithmRecursionSearch(RecursionSearchable<U, V>& generator,
                                 AbstractGenerator<V>& bg) {
            rand = &generator;
            baseGenerator = &bg;
            count = 0;
            isPrime = &MersennePrimitivity;
        }

        /**
         *\japanese
         * コンストラクタ
         *
         * このコンストラクタでは、原始多項式判定アルゴリズムを指定できる。
         *
         * bg はgenerator とは異なるインスタンスである必要がある。
         * 通常の場合、bg にはMersenneTwisterのインスタンスを与えるとよい。
         * TinyMT では、規則的に減少する数列を使用している。
         *
         * @param generator generator はこのクラスのメソッドによって変更される。
         * @param bg パラメータ探索に使用する疑似乱数生成器
         * @param primitivity 原始多項式判定アルゴリズム
         *\endjapanese
         *
         *\english
         * Constructor
         * Users can specify algorithm of primitivity judgement.
         * @param[in,out] generator GF(2)-linear generator whose parameters will
         * be searched. The generator will be changed.
         * @param[in,out] bg a generator used for generating numbers to make
         * parameters. This generator is not need to be GF(2)-linear
         * pseudo random number generator, for example, TinyMTDC in
         * sample directory uses sequential counter.
         * @param[in] primitivity A class to judge primitivity.
         *\endenglish
         */
        AlgorithmRecursionSearch(RecursionSearchable<U, V>& generator,
                                 AbstractGenerator<V>& bg,
                                 const AlgorithmPrimitivity& primitivity) {
            rand = &generator;
            baseGenerator = &bg;
            count = 0;
            isPrime = &primitivity;
        }

        /**
         *\japanese
         * 状態遷移パラメータの探索を開始する
         *<ol>
         * <li>疑似乱数生成器に状態遷移パラメータをランダムに生成させ、
         * <li>その状態遷移パラメータのもとで出力列の最小多項式を求める。
         * <li>最小多項式が求める次数の原始多項式か判定し、原始多項式なら成功して終了。
         * <li>そうでなければ、状態遷移パラメータのランダム生成から繰り返す。
         * <li>繰り返しの回数が try_count を越えると失敗して終了する。
         *</ol>
         * @param try_count 試行回数の上限
         * @return true 求める次数の原始多項式となる最小多項式が得られた場合
         *\endjapanese
         *
         *\english
         * Start searching recursion parameters.
         *<ol>
         * <li>the generator make parameters randomly
         * (RecursionSearchable.setUpParam()),
         * <li>calculate minimal polynomial of the generator,
         * <li>check if the polynomial has max possible degree and
         * is primitive.
         * <li>if check OK, then return true.
         * <li>else repeat from 1.
         *</ol>
         * If repeat count is over \b try_count, then return false.
         * @param[in] try_count maximum count of try
         * @return true when proper minimal polynomial is gotten.
         *\endenglish
         */
        bool start(int try_count) {
            long size = rand->bitSize();
            long degree;
            for (int i = 0; i < try_count; i++) {
                rand->setUpParam(*baseGenerator);
                rand->seed(getOne<U>());
                minpoly(poly, *rand);
                count++;
                degree = deg(poly);
                if (degree != size) {
#if defined(DEBUG)
                    cout << "degree:" << degree << endl;
#endif
                    continue;
                }
                if ((*isPrime)(size, poly)) {
                    return true;
                }
            }
            return false;
        }

        /**
         *\japanese
         * 疑似乱数生成器のパラメータを表す文字列を返す
         * このメソッドは start() が true を返した場合にのみ呼び出すべきである。
         * @return 疑似乱数生成器のパラメータを表す文字列
         *\endjapanese
         *
         *\english
         * Returns a string which shows parameters of pseudo random
         * Call this method only after start() returns true.
         * @return String which shows parameters of pseudo random
         * number generator.
         *\endenglish
         */
        const std::string getParamString() {
            return rand->getParamString();
        }

        /**
         *\japanese
         * 最小多項式を返す
         * このメソッドは start() が true を返した場合にのみ呼び出すべきである。
         *
         * 一般には最小多項式は初期状態と出力関数に依存するが、GF(2)線
         * 形疑似乱数生成器の最大周期を与えるような最小多項式は一つであ
         * る。
         *
         * @return 最小多項式
         *\endjapanese
         *
         *\english
         * Returns a minimal polynomial of output of the pseudo
         * random number generator.
         * Call this method only after start() returns true.
         * @return a minimal polynomial
         *\endenglish
         */
        const NTL::GF2X& getMinPoly() const {
            return poly;
        }

        /**
         *\japanese
         * このインスタンスが作られてからstart() が終了するまでに試行した回数を返す。
         * @return tried count after this class has created.
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
        RecursionSearchable<U, V> *rand;
        AbstractGenerator<V> *baseGenerator;
        const AlgorithmPrimitivity *isPrime;
        NTL::GF2X poly;
        long count;
    };
}
#endif // MTTOOLBOX_ALGORITHM_RECURSION_SEARCH_HPP
