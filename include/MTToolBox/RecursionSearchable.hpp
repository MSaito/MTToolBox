#ifndef MTTOOLBOX_RECURSION_SEARCHABLE_HPP
#define MTTOOLBOX_RECURSION_SEARCHABLE_HPP

/**
 * @file RecursionSearchable.hpp
 *\japanese
 * @brief このクラスは状態遷移パラメータ探索を行えるGF(2)線形疑似乱数生
 * 成器のクラスである。
 *\endjapanese
 *
 *\english
 * @brief Abstract class for searching parameters of state transition
 * function of pseudo random number generator.
 *\endenglish
 *
 *
 */
#include <stdint.h>
#include <inttypes.h>
#include <string>
#include <MTToolBox/AbstractGenerator.hpp>

namespace MTToolBox {
    /**
     * @class ResursionSearchable
     *\japanese
     * @brief このクラスは状態遷移パラメータ探索を行えるGF(2)線形疑似乱
     * 数生成器のクラスである。
     *
     * @tparam U 疑似乱数生成器の出力のタイプ、符号なし型であること
     *\endjapanese
     *
     *\english
     * @brief Abstract class for searching parameters of state transition
     * function of pseudo random number generator.
     * @tparam U type of output of pseudo random number generator, should
     * be unsigned integer.
     *\endenglish
     */
    template<class U>
    class RecursionSearchable : public AbstractGenerator<U> {
    public:

        /**
         *\japanese
         * 仮想デストラクタ（必須）
         *\endjapanese
         *
         *\english
         * Virtual destructor (always required)
         *\endenglish
         */
        virtual ~RecursionSearchable() {};

        /**
         *\japanese
         * 状態遷移関数のパラメータを設定する。
         *
         * このメソッドの内部で、別の疑似乱数生成器を使用してパラメータを
         * 設定すること。TinyMT では疑似乱数ではなくカウントダウンする
         * 数値を元にパラメータを設定しているが、そうしてもよい。
         * @param[in,out] generator 疑似乱数生成器
         *\endjapanese
         *
         *\english
         * Users should set parameters for their generator when
         * this method is called.
         * @param[in, out] generator Source of random parameter set.
         * \b generator may be Mersenne Twister or SequentialGenerator.
         *\endenglish
         */
        virtual void setUpParam(AbstractGenerator<U>& generator) = 0;

        /**
         *\japanese
         * パラメータのヘッダ文字列を返す。
         *
         * パラメータを出力する際に、わかりやすいようにヘッダを表示する。
         * このメソッドは、Dynamic Creator のように大量のパラメータを探
         * 索する際に使用する。ここで返却する文字列パラメータは状態遷移関数の
         * パラメータだけでなく、テンパリングパラメータも含めてよい。
         * @return パラメータのヘッダ文字列
         *\endjapanese
         *
         *\english
         * Returns header string of parameters.
         * @return header string of parameters.
         *\endenglish
         */
        virtual const std::string getHeaderString() = 0;

        /**
         *\japanese
         * パラメータの文字列表現を返す。
         *
         * Dynamic Creator のように大量のパラメータを探索する場合は、
         * getHeaderString() と組み合わせて出力するとよい。そうでない場合は、
         * getHeaderString() は何もしないようにしてもよい。ここで出力するパ
         * ラメータは状態遷移関数のパラメータだけでなく、テンパリングパ
         * ラメータも出力してよい。
         * @return パラメータの文字列
         *\endjapanese
         *
         *\english
         * Returns string expression of parameters.
         * @return string expression of parameters.
         *\endenglish
         */
        virtual const std::string getParamString() = 0;
    };
}

#endif // MTTOOLBOX_RECURSION_SEARCHABLE_HPP
