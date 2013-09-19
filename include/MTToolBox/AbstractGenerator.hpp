#ifndef MTTOOLBOX_ABSTRACT_GENERATOR_HPP
#define MTTOOLBOX_ABSTRACT_GENERATOR_HPP
/**
 * @file AbstractGenerator.hpp
 *
 * @brief GF(2)線形疑似乱数生成器の抽象クラス
 *
 * @author Mutsuo Saito (Hiroshima University)
 * @author Makoto Matsumoto (The University of Tokyo)
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
 * MTToolBox の名前空間
 */
namespace MTToolBox {
    /**
     * @class AbstractGenerator
     * GF(2)線形の状態遷移をもつ疑似乱数生成器
     * @tparam U 疑似乱数生成器の出力の型
     */
    template<class U>
    class AbstractGenerator {
    public:

        /**
         * 内部状態を次状態に遷移し、疑似乱数をひとつ出力する。
         * @return 疑似乱数
         */
        virtual U generate() = 0;

        /**
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
         */
        virtual void seed(U value) = 0;

        /**
         * 内部状態空間のビットサイズを返す。
         *
         * 内部状態空間のビットサイズ、つまりGF(2)線形空間の次元を返す。
         * メルセンヌツイスタのように欠けた配列を使用している場合は、実
         * 際に使用されている部分のビットサイズを返すこと。
         */
        virtual int bitSize() const = 0;
    };
}
#endif //MTTOOLBOX_ABSTRACT_GENERATOR_HPP
