#ifndef MTTOOLBOX_TEMPERING_CALCULATABLE_HPP
#define MTTOOLBOX_TEMPERING_CALCULATABLE_HPP
/**
 * @file TemperingCalculatable.hpp
 * @brief テンパリングパラメータ探索用の抽象クラス
 *
 * テンパリングを行うGF(2)線形疑似乱数生成器は、このクラスを継承する
 * ことによって、TemperingAlgorithmを使用したテンパリングパラメータ
 * 探索が可能になる。
 */
#include <stdint.h>
#include <inttypes.h>
#include <MTToolBox/EquidistributionCalculatable.hpp>

namespace MTToolBox {
    /**
     * @class TemperingCalculatable
     *
     * テンパリングを行うGF(2)線形疑似乱数生成器は、このクラスを継承す
     * ることによって、TemperingAlgorithmを使用したテンパリングパラメー
     * タ探索が可能になる。
     */
    template<class U>
    class TemperingCalculatable : public EquidistributionCalculatable<U> {
    public:

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
        virtual void setUpTempering() = 0;

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
