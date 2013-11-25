#ifndef MTTOOLBOX_SEQUENTIAL_HPP
#define MTTOOLBOX_SEQUENTIAL_HPP
/**
 * @file Sequential.hpp
 *
 *\japanese
 * @brief 順番にカウントダウンする数を生成する
 *\endjapanese
 *
 *\english
 * @brief Generates counting down numbers
 *\endenglish
 *
 * @author Mutsuo Saito (Hiroshima University)
 * @author Makoto Matsumoto (Hiroshima University)
 *
 * Copyright (c) 2013 Mutsuo Saito, Makoto Matsumoto and Hiroshima
 * University. All rights reserved.
 *
 * The 3-clause BSD License is applied to this software, see
 * LICENSE.txt
 */

#include <stdexcept>
#include <stdint.h>
#include <inttypes.h>
#include <MTToolBox/AbstractGenerator.hpp>
#include <MTToolBox/util.hpp>

namespace MTToolBox {
    /**
     * @class Sequential
     *\japanese
     * @brief カウントダウン生成器
     *
     *
     * @tparam T 出力のタイプ、整数型であること
     *\endjapanese
     *
     *\english
     * @brief Counting down number generator
     * @tparam T type of output, should be integer type
     *\endenglish
     */
    template<typename T>
    class Sequential : public AbstractGenerator<T> {
    public:

        /**
         *\japanese
         * コンストラクタ
         *\endjapanese
         *
         *\english
         * Constructor
         *\endenglish
         */
        Sequential() {
        }

        /**
         *\japanese
         * マスク付きコンストラクタ
         * @param[in] p_mask 出力の際にカウンタと排他的論理和を取るためのマスク
         *\endjapanese
         *
         *\english
         * Constructor with mask
         * @param[in] p_mask \b p_mask and internal counter are
         * exclusively or-ed when output time.
         *\endenglish
         */
        Sequential(T p_mask) {
            status = static_cast<T>(-1);
            mask = p_mask;
            error = false;
        }

        /*
         *\japanese
         * マスクとシード付きコンストラクタ
         * @param[in] p_mask 出力の際にカウンタと排他的論理和を取るためのマスク
         * @param[in] seed 内部カウンタの初期値
         *\endjapanese
         *
         *\english
         * Constructor with mask and seed
         * @param[in] seed initial value of the internal counter
         * @param[in] p_mask \b p_mask and internal counter are
         * exclusively or-ed when output time.
         *\endenglish
         */
        Sequential(T p_mask, T seed) {
            status = seed;
            mask = p_mask;
            error = false;
        }

        /*
         *\japanese
         * コピーコンストラクタ
         * @param[in] src コピー元
         *\endjapanese
         *
         *\english
         * Copy Constructor
         * @param[in] src source of copy
         *\endenglish
         */
        Sequential(Sequential<T>& src) : AbstractGenerator<T>() {
            status = src.status;
            mask = src.mask;
            error = src.error;
        }

        /*
         *\japanese
         * 初期化
         * @param[in] value 内部カウンタの初期値
         *\endjapanese
         *
         *\english
         * Initialization
         * @param[in] value Initial value of the internal counter
         *\endenglish
         */
        void seed(T value) {
            reseed(value);
        }

        /*
         *\japanese
         * 初期化
         * @param[in] seed 内部カウンタの初期値
         *\endjapanese
         *
         *\english
         * Initialization
         * @param[in] seed Initial value of the internal counter
         *\endenglish
         */
        void reseed(T seed) {
            status = seed;
            error = false;
        }

        /*
         *\japanese
         * 次の数を返す
         *
         * このメソッドはnext()を呼び出している。
         * @see next()
         * @return next value
         *\endjapanese
         *
         *\english
         * Returns next value
         *
         * This method calls next()
         * @see next()
         * @return next value
         *\endenglish
         */
        T generate() {
            return next();
        }

        /*
         *\japanese
         * 次の数を返す
         *
         * 内部カウンタとマスクとの排他的論理和をとって返す。
         * 返却値を決定後に、内部カウンタをひとつ減らす。
         * @throw std::underflow_exception ゼロを返した後、さらにこのメ
         * ソッドが呼ばれた場合
         *\endjapanese
         *
         *\english
         * Returns next value
         *
         * Returns exclusive or of the internal counter and
         * mask given by constructor argument.
         * After return value is decided, the internal counter
         * will be decremented.
         * @throws std::underflow_exception when this method
         * is called after this method returns zero.
         * @return next value
         *\endenglish
         */
        T next() {
            if (error) {
                throw std::underflow_error("count over zero exception");
            }
            if (status <= 0) {
                error = true;
            }
            T work = status;
            status -= 1;
            return work ^ mask;
        }

        /*
         *\japanese
         * 内部カウンタのビットサイズを返す。
         *
         * AbstracutGenerator に合わせるため
         * @return 内部カウンタのビットサイズ
         *\endjapanese
         *
         *\english
         * Returns bit size of the internal counter
         *
         * To fit to interface of AbstructGenerator
         * @return bit size of internal counter
         *\endenglish
         */
        int bitSize() const {
            int r = bit_size<T>();
            return r;
        }
    private:
        T status;
        T mask;
        bool error;
    };
}

#endif //  MTTOOLBOX_SEQUENTIAL_HPP


