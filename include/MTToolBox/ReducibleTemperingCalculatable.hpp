#ifndef MTTOOLBOX_REDUCIBLE_TEMPERING_CALCULATABLE_HPP
#define MTTOOLBOX_REDUCIBLE_TEMPERING_CALCULATABLE_HPP
/**
 * @file ReducibleTemperingCalculatable.hpp
 *\japanese
 * @brief 可約ジェネレータのテンパリングパラメータ探索用の抽象クラス
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
 * @author Mutsuo Saito (Manieth Corp.)
 * @author Makoto Matsumoto (Hiroshima University)
 *
 * Copyright (C) 2015 Mutsuo Saito, Makoto Matsumoto, Manieth Corp.
 * and Hiroshima University.
 * All rights reserved.
 *
 * The 3-clause BSD License is applied to this software, see
 * LICENSE.txt
 */
#include <stdint.h>
#include <inttypes.h>
#include <MTToolBox/ReducibleGenerator.hpp>
#include <MTToolBox/TemperingCalculatable.hpp>

namespace MTToolBox {
    /**
     * @class TemperingCalculatable
     *
     *\japanese
     * テンパリングを行う可約ジェネレータは、このクラスを継承す
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
    class ReducibleTemperingCalculatable
        : virtual public ReducibleGenerator<U>,
          virtual public TemperingCalculatable<U> {
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
        virtual ~ReducibleTemperingCalculatable() {}
    };
}

#endif // MTTOOLBOX_TEMPERING_CALCULATABLE_HPP
