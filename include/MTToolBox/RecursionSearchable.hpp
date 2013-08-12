#ifndef MTTOOLBOX_RECURSION_SEARCHABLE_HPP
#define MTTOOLBOX_RECURSION_SEARCHABLE_HPP

/**
 * @file RecursionSearchable.hpp
 * @brief このクラスは状態遷移パラメータ探索を行えるGF(2)線形疑似乱数生成器のクラスである。
 *
 *
 */
#include <stdint.h>
#include <inttypes.h>
#include <MTToolBox/AbstractGenerator.hpp>

namespace MTToolBox {
    /**
     * @class ResursionSearchable
     * @brief このクラスは状態遷移パラメータ探索を行える
     * GF(2)線形疑似乱数生成器のクラスである。
     *
     * @tparam U 疑似乱数生成器の出力のタイプ
     */
    template<class U>
    class RecursionSearchable : public AbstractGenerator<U> {
    public:

        /**
         * 状態遷移関数のパラメータを設定する。
         *
         * このメソッドの内部で、別の疑似乱数生成器を使用してパラメータを
         * 設定すること。TinyMT では疑似乱数ではなくカウントダウンする
         * 数値を元にパラメータを設定しているが、そうしてもよい。
         */
        virtual void setUpParam() = 0;

        /**
         * パラメータのヘッダを出力する。
         *
         * パラメータを出力する際に、わかりやすいようにヘッダを表示する。
         * このメソッドは、Dynamic Creator のように大量のパラメータを探
         * 索する際に使用する。ここで出力するパラメータは状態遷移関数の
         * パラメータだけでなく、テンパリングパラメータも出力してよい。
         */
        virtual void printHeader(std::ostream& out) = 0;

        /**
         * パラメータを出力する。
         *
         * Dynamic Creator のように大量のパラメータを探索する場合は、
         * printHeader() と組み合わせて出力するとよい。そうでない場合は、
         * printHeader() は何もしないようにしてもよい。ここで出力するパ
         * ラメータは状態遷移関数のパラメータだけでなく、テンパリングパ
         * ラメータも出力してよい。
         */
        virtual void printParam(std::ostream& out) = 0;
    };
}

#endif // MTTOOLBOX_RECURSION_SEARCHABLE_HPP
