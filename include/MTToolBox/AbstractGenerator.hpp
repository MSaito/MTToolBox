#ifndef MTTOOLBOX_ABSTRACT_GENERATOR_HPP
#define MTTOOLBOX_ABSTRACT_GENERATOR_HPP
/**
 * @file AbstractGenerator.hpp
 *
 *\japanese
 * @brief GF(2)線形疑似乱数生成器の抽象クラス
 *\endjapanese
 *\english
 * @brief Abstract class of GF(2)-linear pseudo random number generators.
 *\endenglish
 *
 * @author Mutsuo Saito (Hiroshima University)
 * @author Makoto Matsumoto (Hiroshima University)
 *
 * Copyright (C) 2013 Mutsuo Saito, Makoto Matsumoto
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
     * @class AbstractGenerator
     *\japanese
     * 疑似乱数生成器
     * @tparam U 疑似乱数生成器の出力の型
     *\endjapanese
     *
     *\english
     * pseudo random number generator.
     * @tparam U output type of a pseudo random number generator.
     *\endenglish
     */
    template<class U>
    class AbstractGenerator {
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
        virtual ~AbstractGenerator() {}

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
        virtual U generate() = 0;

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
        virtual void seed(U value) = 0;

        /**
         *\japanese
         * 内部状態空間のビットサイズを返す。
         *
         * 内部状態空間のビットサイズ、つまりGF(2)線形空間の次元を返す。
         * メルセンヌツイスタのように欠けた配列を使用している場合は、実
         * 際に使用されている部分のビットサイズを返すこと。
         * @return 状態空間のビットサイズ
         *\endjapanese
         *
         *\english
         *
         * Return bit size of internal state, i.e dimension of
         * GF(2)-vector space. It will be Mersenne Exponent, when
         * generator use incomplete array, like Mersenne Twister.
         *
         * @return bit size of internal state
         *\endenglish
         */
        virtual int bitSize() const = 0;
    };
}
#endif //MTTOOLBOX_ABSTRACT_GENERATOR_HPP
