#ifndef MTTOOLBOX_ALGORITHM_RECURSION_SEARCH_HPP
#define MTTOOLBOX_ALGORITHM_RECURSION_SEARCH_HPP
/**
 * @file AlgorithmRecursionSearch.hpp
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
#include <MTToolBox/AlgorithmPrimitivity.hpp>
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
         * このコンストラクタでは、状態空間の大きさがメルセンヌ指数の場合を
         * 想定している。
         *
         * bg はgenerator とは異なるインスタンスである必要がある。
         * 通常の場合、bg にはMersenneTwisterのインスタンスを与えるとよい。
         * TinyMT では、規則的に減少する数列を使用している。
         *
         * @param generator generator はこのクラスのメソッドによって変更される。
         * @param bg パラメータ探索に使用する疑似乱数生成器
         */
        AlgorithmRecursionSearch(RecursionSearchable<U>& generator,
                                 AbstractGenerator<U>& bg) {
            rand = &generator;
            baseGenerator = &bg;
            count = 0;
            isPrime = &MersennePrimitivity;
        }

        /**
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
         */
        AlgorithmRecursionSearch(RecursionSearchable<U>& generator,
                                 AbstractGenerator<U>& bg,
                                 const AlgorithmPrimitivity& primitivity) {
            rand = &generator;
            baseGenerator = &bg;
            count = 0;
            isPrime = &primitivity;
        }

        /**
         * 疑似乱数生成器に状態遷移パラメータをランダムに生成させ、
         * その状態遷移パラメータのもとで出力列の最小多項式を求める。
         * 最小多項式が求める次数の原始多項式か判定し、原始多項式なら成功して終了。
         * そうでなければ、状態遷移パラメータのランダム生成から繰り返す。
         * 繰り返しの回数が try_count を越えると失敗して終了する。
         *
         * @param try_count 試行回数の上限
         * @return true 求める次数の原始多項式となる最小多項式が得られた場合
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
                if ((*isPrime)(size, poly)) {
                    return true;
                }
            }
            return false;
        }

        /**
         * このメソッドは start() が true を返した場合にのみ呼び出すべきである。
         * @return 疑似乱数生成器のパラメータを表す文字列
         */
        const std::string getParamString() {
            return rand->getParamString();
        }

        /**
         * このメソッドは start() が true を返した場合にのみ呼び出すべきである。
         *
         * 一般には最小多項式は初期状態と出力関数に依存するが、GF(2)線
         * 形疑似乱数生成器の最大周期を与えるような最小多項式は一つであ
         * る。
         *
         * @return 最小多項式
         */
        const NTL::GF2X& getMinPoly() const {
            return poly;
        }

        /**
         * このインスタンスが作られてからstart() が終了するまでに試行した回数を返す。
         * @return tried count after this class has created.
         */
        long getCount() const {
            return count;
        }

    private:
        RecursionSearchable<U> *rand;
        AbstractGenerator<U> *baseGenerator;
        const AlgorithmPrimitivity *isPrime;
        NTL::GF2X poly;
        long count;
    };
}
#endif // MTTOOLBOX_ALGORITHM_RECURSION_SEARCH_HPP
