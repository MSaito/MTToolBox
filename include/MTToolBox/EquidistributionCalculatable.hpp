#ifndef MTTOOLBOX_EQUIDISTRIBUTION_CALCULATABLE_HPP
#define MTTOOLBOX_EQUIDISTRIBUTION_CALCULATABLE_HPP
/**
 * @file EquidistributionCalculatable.hpp
 * @brief このクラスはGF(2)線形疑似乱数生成器の均等分布次元を計算するためのクラスである。
 *
 *
 */
#include <stdint.h>
#include <inttypes.h>
#include <MTToolBox/RecursionSearchable.hpp>

namespace MTToolBox {
    /**
     * @class EquidistributionCalculatable
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
     * @tparam U 疑似乱数生成器の出力のタイプ
     */
    template<class U>
    class EquidistributionCalculatable : public RecursionSearchable<U> {
    public:
#if 0
        /**
         * 仮想デストラクタ
         */
        virtual ~EquidistributionCalculatable<U>() = 0;
#endif
        /*
         * 自分のコピーを返す。
         *
         * @note clone には問題があるが、テンプレートよりよいと判断する。
         * @return 自分自身のコピー
         */
        virtual EquidistributionCalculatable<U> * clone() const = 0;

        /*
         * 上位(MSBから) @code outBitLen だけ出力する。
         *
         * それ以外のビットは0でなければならない。
         * 均等分布次元の計算にとって重要なメソッドである。
         * このメソッドを実行すると、状態空間の遷移も行われる。
         * @param[in] outBitLen 上位から何ビット出力するかを指定する。
         * @return 上位 @code outBitLen 以外は0の出力
         */
        virtual U generate(int outBitLen) = 0;

        /**
         * GF(2)線形疑似乱数生成器の状態空間を加算し、内部状態を変更する。
         *
         * @caution 配列を使ってラウンドロビン方式で状態空間を保持している場合、
         * このメソッドを実装する際に、インデックスを一致させて加算することが重要である。
         * \b state を状態配列、\b size を配列のサイズ、\b index をラウンドロビンの
         * インデックスとすると以下のような加算が必要になる。
         * @verbatim
           for (int i = 0; i < size; i++) {
             state[(index + i) % size] ^= that.state[(that.index + i) % size];
           }
          @endverbatim
         * これはもう少し効率的に書くことが出来る。（サンプルを参照）
         */
        virtual void add(EquidistributionCalculatable& that) = 0;
#if 0
        /**
         * 状態空間を次状態に遷移する
         *
         */
        virtual void next() = 0;
#endif
        /**
         * 状態空間をすべてゼロにセットする。
         */
        virtual void setZero() = 0;

        /**
         * 状態空間がすべてゼロかチェックする。
         * @caution メルセンヌツイスタのように、欠けた配列を使用する場合は、
         * 真の状態空間がすべてゼロかチェックする必要がある。言い換えると、
         * 必要なマスクをして余分な部分を排除してからゼロチェックする
         * べきである。
         *
         * @return true なら状態空間がすべてゼロ
         */
        virtual bool isZero() = 0;
#if 0
        virtual void setUpParam() = 0;
        virtual void printParam(std::ostream& out) = 0;
#endif
    };
}

#endif // MTTOOLBOX_EQUIDISTRIBUTION_CALCULATABLE_HPP
