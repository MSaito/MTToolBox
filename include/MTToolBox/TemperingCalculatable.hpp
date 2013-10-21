#ifndef MTTOOLBOX_TEMPERING_CALCULATABLE_HPP
#define MTTOOLBOX_TEMPERING_CALCULATABLE_HPP
/**
 * @file TemperingCalculatable.hpp
 *\japanese
 * @brief テンパリングパラメータ探索用の抽象クラス
 *
 * テンパリングを行うGF(2)線形疑似乱数生成器は、このクラスを継承する
 * ことによって、TemperingAlgorithmを使用したテンパリングパラメータ
 * 探索が可能になる。
 *\endjapanese
 *
 *\english
 * @brief Abstruct class for searching tempering parameters.
 *
 * Users can search tempering parameters by making GF(2)-linear
 * pseudo random generator class which inherits from this class.
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
#include <stdint.h>
#include <inttypes.h>
#include <MTToolBox/EquidistributionCalculatable.hpp>

namespace MTToolBox {
    /**
     * @class TemperingCalculatable
     *
     *\japanese
     * テンパリングを行うGF(2)線形疑似乱数生成器は、このクラスを継承す
     * ることによって、TemperingAlgorithmを使用したテンパリングパラメー
     * タ探索が可能になる。
     * @tparam U 疑似乱数生成器の出力の型、符号なし型であること
     *\endjapanese
     *
     *\english
     * Users can search tempering parameters by making GF(2)-linear
     * pseudo random generator class which inherits from this class.
     * @tparam U type of output of pseudo random number generator
     *\endenglish
     */
    template<class U>
    class TemperingCalculatable : public EquidistributionCalculatable<U> {
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
        virtual ~TemperingCalculatable() {}

        /**
         *\japanese
         * テンパリングパラメータをセットする。
         * @param[in] mask \b pattern のうち実際にテンパリングパラメータに
         * 設定するべき部分を指定する。
         * @param[in] pattern テンパリングパラメータの一部にセットされるパターン
         * @param[in] index テンパリングパラメータの数が1より大きいときに
         * 何番目のテンパリングパラメータかを示す。
         *\endjapanese
         *
         *\english
         * Set tempering parameters
         * @param[in] mask Shows which part of \b pattern should be
         * set to parameter.
         * @param[in] pattern a part of pattern set to tempering parameter
         * @param[in] index index of tempering parameters.
         *\endenglish
         */
        virtual void setTemperingPattern(U mask, U pattern, int index) = 0;

        /**
         *\japanese
         * テンパリングテーブルの準備が必要な場合はここで準備する。
         * MTGP の場合はルックアップテーブルの準備をしている。
         * テンパリングパラメータの数が1の場合はおそらく準備をする必要はない。
         *\endjapanese
         *
         *\english
         * If preparing is needed before generation, here is the place
         * to prepare tempering parameters.
         *\endenglish
         */
        virtual void setUpTempering() = 0;

        /**
         *\japanese
         * 出力の上位ビットと下位ビットを反転する。
         * これは下位ビットから見た均等分布次元を計測またはテンパリング
         * するために使われる。
         *\endjapanese
         *
         *\english
         * Changes bit order of output.
         *
         * MSB of output becomes LSB of output. This is useful
         * to calculate lower bit equi-distribution.
         *\endenglish
         */
        virtual void setReverseOutput() = 0;

        /**
         *\japanese
         * 出力の上位ビットと下位ビットの並びを元に戻す。
         *\endjapanese
         *
         *\english
         * Reset bit order of output.
         *\endenglish
         */
        virtual void resetReverseOutput() = 0;

        /**
         *\japanese
         * 出力の上位ビットと下位ビットが反転しているかを返す。
         * @return true 出力の上位ビットと下位ビットが反転している。
         *\endjapanese
         *
         *\english
         * Shows if bit order is reversed
         * @return true if bit order id reversed
         *\endenglish
         */
        virtual bool isReverseOutput() = 0;
    };
}

#endif // MTTOOLBOX_TEMPERING_CALCULATABLE_HPP
