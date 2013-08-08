#ifndef MTTOOLBOX_TEMPERING_CALCULATABLE_HPP
#define MTTOOLBOX_TEMPERING_CALCULATABLE_HPP
/* -*- coding:utf-8 -*- */
/**
 * @file TemperingCalculatable.hpp
 * @brief テンパリングパラメータ探索用の抽象クラス
 *
 * テンパリングパラメータを探索するアルゴリズム TemperingAlgorithm
 * を利用してテンパリングパラメータを探索する時に使用するクラス
 */
#include <stdint.h>
#include <inttypes.h>
#include <MTToolBox/EquidistributionCalculatable.hpp>

namespace MTToolBox {
    /**
     * @class TemperingCalculatable
     *
     */
    template<class U>
    class TemperingCalculatable : public EquidistributionCalculatable<U> {
    public:
        /**
         * テンパリングパラメータの数を返す。
         * 通常は1を返すようにすればよい。
         * MTGP のようにテーブルを参照してテンパリングする場合は
         * 1より大きな数を返すようにoverrideする。
         * @return テンパリングパラメータの数
         */
        virtual int numberOfTemperingParam() {return 1;}

        /**
         * テンパリングパラメータをセットする。
         * @param[in] mask pattern のうち実際にテンパリングパラメータに
         * 設定するべき部分を指定する。
         * @param[in] pattern テンパリングパラメータの一部にセットされるパターン
         * @param[in] index テンパリングパラメータの数が1より大きいときに
         * 何番目のテンパリングパラメータかを示す。
         */
        virtual void setTemperingPattern(U mask, U pattern, int index) = 0;

        /**
         * テンパリングテーブルの準備が必要な場合はここで準備する。
         * MTGP の場合はルックアップテーブルの準備をしている。
         * テンパリングパラメータの数が1の場合はおそらく準備をする必要はない。
         */
        virtual void setUpTempering() {};

        /**
         * 出力の上位ビットと下位ビットを反転する。
         * これは下位ビットから見た均等分布次元を計測またはテンパリング
         * するために使われる。
         */
        virtual void setReverseOutput() = 0;

        /**
         * 出力の上位ビットと下位ビットの並びを元に戻す。
         */
        virtual void resetReverseOutput() = 0;

        /**
         * 出力の上位ビットと下位ビットが反転しているかを返す。
         * @return true 出力の上位ビットと下位ビットが反転している。
         */
        virtual bool isReverseOutput() = 0;
    };
}

#endif // MTTOOLBOX_TEMPERING_CALCULATABLE_HPP
