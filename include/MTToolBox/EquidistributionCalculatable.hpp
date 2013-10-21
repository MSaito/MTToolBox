#ifndef MTTOOLBOX_EQUIDISTRIBUTION_CALCULATABLE_HPP
#define MTTOOLBOX_EQUIDISTRIBUTION_CALCULATABLE_HPP
/**
 * @file EquidistributionCalculatable.hpp
 *\japanese
 * @brief このクラスはGF(2)線形疑似乱数生成器の均等分布次元を計算するためのクラスである。
 *\endjapanese
 *
 *\english
 * @brief This class is an Abstract class for calculating dimension of
 * equi-distribution for GF(2)-linear pseudo random number generators.
 *\endenglish
 *
 *
 */
#include <stdint.h>
#include <inttypes.h>
#include <MTToolBox/RecursionSearchable.hpp>

namespace MTToolBox {
    /**
     * @class EquidistributionCalculatable
     *\japanese
     * @brief このクラスはGF(2)線形疑似乱数生成器の均等分布次元を計算するためのクラスである。
     *
     * GF(2)線形疑似乱数生成器を設計するならば、均等分布次元を計算して
     * おいた方がよいであろう。均等分布次元を計算するためには、出力関数
     * も含めてGF(2)線形であることが必要である。均等分布次元は状態空間
     * の大きさによって上限が決まるので、状態空間の大きさを固定した時に
     * どのくらい上限に近づけるかという点に理論的な関心がある。一方、多
     * 変数（多次元）の数値積分などの実際の問題領域においては、均等分布
     * 次元の大きさそのものが問題であり、その点から言えば、状態空間の大
     * きな疑似乱数生成器を開発し、均等分布次元の値を大きくすることの方
     * が理論的上限に近いかどうかよりも重要である。
     *
     * @tparam U 疑似乱数生成器の出力のタイプ、符号なし型でなければならない。
     *\endjapanese
     *
     *\english
     * @brief This class is an Abstract class for calculating dimension of
     * equi-distribution for GF(2)-linear pseudo random number generators.
     *
     * @tparam U type of output of pseudo random number generator, should
     * be unsinged number.
     *\endenglish
     */
    template<class U>
    class EquidistributionCalculatable : public RecursionSearchable<U> {
    public:
        using AbstractGenerator<U>::generate;

        /**
         *\japanese
         * 仮想デストラクタ（必須）
         *\endjapanese
         *
         *\english
         * Virtual destructor (always required)
         *\endenglish
         */
        virtual ~EquidistributionCalculatable(){};

        /**
         *\japanese
         * 自分のコピーを返す。
         *
         * @note clone には問題があるが、テンプレートよりよいと判断する。
         * @return 自分自身のコピー
         *\endjapanese
         *
         *\english
         * Return copy of myself.
         * @return copy of myself.
         *\endenglish
         */
        virtual EquidistributionCalculatable<U> * clone() const = 0;

        /**
         *\japanese
         * 上位(MSBから) \b outBitLen だけ出力する。
         *
         * それ以外のビットは0でなければならない。
         * 均等分布次元の計算にとって重要なメソッドである。
         * このメソッドを実行すると、状態空間の遷移も行われる。
         * @param[in] outBitLen 上位から何ビット出力するかを指定する。
         * @return 上位 \b outBitLen 以外は0の出力
         *\endjapanese
         *
         *\english
         * output \b outBitLen from MSB.
         *
         * Remained bits should be all zeros.
         * Users should call statte transition function inside
         * this method.
         * @param[in] outBitLen Shows how many bits are outputed.
         * @return an integer whose \b outBiteLen bits from MSB.
         *\endenglish
         */
        virtual U generate(int outBitLen) = 0;

        /**
         *\japanese
         * GF(2)線形疑似乱数生成器の状態空間を加算し、内部状態を変更する。
         *
         * \warning 配列を使ってラウンドロビン方式で状態空間を保持している場合、
         * このメソッドを実装する際に、インデックスを一致させて加算することが重要である。
         * \b state を状態配列、\b size を配列のサイズ、\b index をラウンドロビンの
         * インデックスとすると以下のような加算が必要になる。
         *\endjapanese
         *
         *\english
         * Add internal state of GF(2)-linear pseudo random number generators.
         *
         * \warning When keeping internal state in an array and using index
         * for round robbin, adding shoud be relative to index like:
         *\endenglish
         * @verbatim
         for (int i = 0; i < size; i++) {
           state[(index + i) % size] ^= that.state[(that.index + i) % size];
         }
         @endverbatim
         */
        virtual void add(EquidistributionCalculatable& that) = 0;

        /**
         *\japanese
         * 状態空間をすべてゼロにセットする。
         *\endjapanese
         *
         *\english
         * Set all zero to internal state.
         *\endenglish
         */
        virtual void setZero() = 0;

        /**
         *\japanese
         * 状態空間がすべてゼロかチェックする。
         *
         * \warning メルセンヌツイスタのように、欠けた配列を使用する場合は、
         * 真の状態空間がすべてゼロかチェックする必要がある。言い換えると、
         * 必要なマスクをして余分な部分を排除してからゼロチェックする
         * べきである。
         *
         * @return true なら状態空間がすべてゼロ
         *\endjapanese
         *
         *\english
         * Check if bits in internal state are all zero.
         *
         * When using an incomplete array like Mersenne Twister,
         * user should check bits in true internal state.
         *\endenglish
         */
        virtual bool isZero() const = 0;
    };
}

#endif // MTTOOLBOX_EQUIDISTRIBUTION_CALCULATABLE_HPP
