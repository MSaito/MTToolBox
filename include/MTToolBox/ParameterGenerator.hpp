#ifndef MTTOOLBOX_PARAMETER_GENERATOR_HPP
#define MTTOOLBOX_PARAMETER_GENERATOR_HPP
/**
 * @file ParameterGenerator.hpp
 *
 *\japanese
 * @brief パラメータ生成器の抽象クラス
 *\endjapanese
 *\english
 * @brief Abstract class of Parameter Generator.
 *\endenglish
 *
 * @author Mutsuo Saito
 * @author Makoto Matsumoto (Hiroshima University)
 *
 * Copyright (C) 2017 Mutsuo Saito, Makoto Matsumoto
 * and Hiroshima University.
 * All rights reserved.
 *
 * The 3-clause BSD License is applied to this software, see
 * LICENSE.txt
 */

/**
 * @namespace MTToolBox
 *
 *\japanese
 * MTToolBox の名前空間
 *\endjapanese
 *\english
 * name space for MTToolBox
 *\endenglish
 */
namespace MTToolBox {
    /**
     * @class ParameterGenerator
     *\japanese
     * 疑似乱数生成器
     *\endjapanese
     *
     *\english
     * pseudo random number generator.
     *\endenglish
     */
    class ParameterGenerator {
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
        virtual ~ParameterGenerator() {}

        /**
         *\japanese
         * 内部状態を次状態に遷移し、疑似乱数をひとつ出力する。
         * @return 疑似乱数
         *\endjapanese
         *
         *\english
         * transit current internal state to next state and output a
         * pseudo random number.
         * @return a pseudo random number
         *\endenglish
         */
        virtual uint32_t getUint32() = 0;

        /**
         *\japanese
         * 内部状態を次状態に遷移し、疑似乱数をひとつ出力する。
         * @return 疑似乱数
         *\endjapanese
         *
         *\english
         * transit current internal state to next state and output a
         * pseudo random number.
         * @return a pseudo random number
         *\endenglish
         */
        virtual uint64_t getUint64() = 0;

        /**
         *\japanese
         * 内部状態を初期化する。
         *
         * 初期化処理は、GF(2)線形である必要はなく、むしろ非GF(2)線形で
         * あることが望ましい。ただし、MTToolBox による周期の計算や均等
         * 分布次元の計算においては、初期化の結果内部状態がゼロでないと
         * いうことだけが重要である。メルセンヌツイスタのように欠けた配
         * 列を使用している場合は、実際に使用している部分がゼロでないよ
         * うに初期化する必要がある。
         *
         * @param[in] value 初期化の種
         *\endjapanese
         *
         *\english
         * initialize internal state
         *
         * Initialization function does not need to be a GF(2)-linear
         * function, non-GF(2)-linear function will be suitable.  But
         * as far as MTToolBox concerned, initialization only need to
         * assure that the internal state is not zero.  When generator
         * use incomplete array, like Mersenne Twister, really used
         * part of array should be set non-zero.
         *
         * @param[in] value seed of initialization
         *\endenglish
         */
        virtual void seed(uint64_t value) = 0;

    };
}
#endif //MTTOOLBOX_PARAMETER_GENERATOR_HPP
