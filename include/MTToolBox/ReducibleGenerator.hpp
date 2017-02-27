#ifndef MTTOOLBOX_REDUCIBLE_GENERATOR_HPP
#define MTTOOLBOX_REDUCIBLE_GENERATOR_HPP
/**
 * @file EquidistributionCalculatable.hpp
 *\japanese
 * @brief このクラスは特性多項式が可約なGF(2)線形疑似乱数生成器を開発するためのクラスである。
 *
 * 状態遷移関数の特性多項式が大きなメルセンヌ指数次の既約因子をもつ
 * GF(2)線形疑似乱数生成器をMTToolBoxでは可約ジェネレータと呼ぶ。
 *
 * 大きな既約因子を持つかどうかの判定や、大きな既約因子に基づく周期を保
 * 証するためのパリティベクターを計算するために必要となるメソッドを備え
 * る。
 *
 * 初めに出力のサイズbsを決める。通常は32ビットまたは64ビットであるが、
 * SIMDを使用する場合はもっと大きくなることもある。基本的には、bs の数
 * 倍より大きなメルセンヌ指数をひとつ選んでMEXPとする。要素のサイズを
 * bs、配列の大きさを MEXP / bs + 1 とする。(MEXP / bs + 1 * bs) が可約
 * ジェネレータの状態空間の大きさとなる。従って状態空間の大きさはbsと
 * MEXPから自動的に決定され、それ以外のサイズは許されない。
 *
 *\endjapanese
 *
 *\english
 * @brief This class is an Abstract class for reducible generator.
 *\endenglish
 *
 *
 */
#include <stdint.h>
#include <inttypes.h>
#include <NTL/GF2X.h>
#include <MTToolBox/EquidistributionCalculatable.hpp>

namespace MTToolBox {
    /**
     * @class ReducibleGenerator
     *\japanese
     * @brief このクラスは特性多項式が可約なGF(2)線形疑似乱数生成器を開
     * 発するためのクラスである。
     *
     * 大きな既約因子を持つかどうかの判定や、大きな既約因子に基づく周期を保
     * 証するためのパリティベクターを計算するために必要となるメソッドを備え
     * る。
     * @tparam U 疑似乱数生成器の出力の型
     *\endjapanese
     *
     *\english
     * @brief This class is an Abstract class for reducible generator.
     * @tparam U output type of the generator.
     *\endenglish
     *
     */
    template<typename U>
    class ReducibleGenerator
        : virtual public EquidistributionCalculatable<U> {
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
        virtual ~ReducibleGenerator(){};

        /**
         *\japanese
         * 状態空間の中の指定された1ビットだけ1にし、それ以外をすべてゼロにする。
         * このメソッドはパリティベクトルを計算するために使用される。
         *@param[in] bitPos 1 にセットするビット位置
         *\endjapanese
         *
         *\english
         * set 1 for specified bit and clear all other bits in internal state.
         * This method is used for calculating parity vector.
         *@param[in] bitPos position of bit set to be one
         *\endenglish
         */
        virtual void setOneBit(int bitPos) = 0;

        /**
         *\japanese
         * パリティチェック用の状態空間の位置にある値を返す。
         * @return パリティチェック用の位置に今ある値
         *\endjapanese
         *
         *\english
         * Virtual destructor (always required)
         * Returs a value in the position for parity check
         * @return a value in the position for parity check
         *\endenglish
         */
        virtual U getParityValue() const = 0;

        /**
         *\japanese
         * パリティチェック用の状態空間の位置に指定された値をセットする。
         * @param[in] parity セットする値
         *\endjapanese
         *
         *\english
         * Set a specified value in the position for parity check.
         * @param[in] parity a value to be set
         *\endenglish
         */
        virtual void setParityValue(U parity) = 0;

        /**
         *\japanese
         * この疑似乱数生成器の最低周期となるメルセンヌ指数を取得する。
         * @return メルセンヌ指数
         *\endjapanese
         *
         *\english
         * Get Mersenne Exponent which is a exponent part of certified
         * minimum period.
         * @return mersenne exponent
         *\endenglish
         */
        virtual int getMexp() const = 0;

    };

    /**
     *\japanese
     * 可約疑似乱数生成器の状態空間を多項式で殲滅する。
     * @tparam U 疑似乱数生成器の出力の型
     * @param[in,out] rg 可約疑似乱数生成器
     * @param[in] poly 殲滅多項式
     *\endjapanese
     *
     *\english
     * Annihilate internal state of generator by given polynomial.
     * @tparam U output type of the generator.
     * @param[in,out] rg reducible generator
     * @param[in] poly annihilator polynomial
     *\endenglish
     */
    template<typename U>
    void annihilate(EquidistributionCalculatable<U>* rg,
                    const NTL::GF2X& poly) {
        EquidistributionCalculatable<U> *other = rg->clone();
        rg->setZero();
        for (int i = 0; i <= deg(poly); i++) {
            if (coeff(poly, i) != 0) {
                rg->add(*other);
            }
            other->generate();
        }
        delete other;
    }
}

#endif // MTTOOLBOX_REDUCIBLE_GENERATOR_HPP
