#ifndef MTTOOLBOX_ALGORITHM_EQUIDISTRIBUTION_HPP
#define MTTOOLBOX_ALGORITHM_EQUIDISTRIBUTION_HPP
/**
 * @file AlgorithmEquidistribution.hpp
 *
 * @brief 疑似乱数生成器の出力の均等分布次元を計算する。
 *
 * PIS法[1](原瀬)によって疑似乱数生成器の出力の均等分布次元を計算するアルゴリズム
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
 *
 * [1] S. Harase. An efficient lattice reduction method for F2-linear
 * pseudorandom number generators using Mulders and Storjohann
 * algorithm. Journal of Computational and Applied Mathematics,
 * 236(2):141–149, August 2011. doi:10.1016/j.cam.2011.06.005.
 */
#include <tr1/memory>
#include <stdexcept>
#include <MTToolBox/EquidistributionCalculatable.hpp>
#include <MTToolBox/util.hpp>

namespace MTToolBox {
    /**
     * @class linear_generator_vector
     * @brief GF(2)ベクトルとしてのGF(2)疑似乱数生成器
     *
     * このクラスはGF(2)ベクトルとしての疑似乱数生成器を表す。また
     * GF(2)係数多項式をv個まとめたベクトルとしても扱うが、多項式特有の
     * 処理としては不定元倍(x をかける）ことのみをサポートする。これは
     * next_state()メソッドによって行われる。
     *
     * @tparam T 疑似乱数生成器の出力の型
     */
    template<typename T>
    class linear_generator_vector {
    public:

        /**
         * 均等分布次元計算可能な疑似乱数生成器
         */
        typedef EquidistributionCalculatable<T> ECGenerator;

        /**
         * コンストラクタ
         *
         * @param generator GF(2)疑似乱数生成器
         */
        linear_generator_vector<T>(const ECGenerator& generator) {
            using namespace std::tr1;
            shared_ptr<ECGenerator>
                r(reinterpret_cast<ECGenerator *>(generator.clone()));
            rand = r;
            rand->seed(1);
            count = 0;
            zero = false;
            next = 0;
        }

        /**
         * 標準基底コンストラクタ
         *
         * 標準基底を構成するベクトルのひとつを生成する。
         * 出力の特定のビットに1度だけ1を出力して、その後はずっと0を出力する
         * 疑似乱数生成器のコンストラクタ。例えば、1, 0, 0, 0, .... と出力する、
         * あるいは 8, 0, 0, 0, .... と出力するなど。
         *
         * @param generator GF(2)疑似乱数生成器
         * @param bit_pos 1 の位置
         */
        linear_generator_vector<T>(const ECGenerator& generator,
                                   int bit_pos) {
            using namespace std::tr1;
            shared_ptr<ECGenerator>
                r(reinterpret_cast<ECGenerator *>(generator.clone()));
            rand = r;
            rand->setZero();
            count = 0;
            zero = false;
            next = static_cast<T>(1) << (sizeof(T) * 8 - bit_pos - 1);
        }

        void add(const linear_generator_vector<T>& src);
        void next_state(int bit_len);
        void debug_print();

        /**
         * GF(2)線形疑似乱数生成器
         */
        std::tr1::shared_ptr<ECGenerator> rand;

        /**
         * next_state() が呼ばれた回数
         * これは多項式としてみた場合の次数に関係がある。
         */
        int count;

        /**
         * ゼロベクトルであるかどうかを示す。
         */
        bool zero;

        /**
         * 疑似乱数生成器の最新の出力（の上位vビット）または
         * 多項式の最高次の係数のなすベクトル
         */
        T next;
    };

    /**
     * @class AlgorithmEquidistribution
     * @brief 疑似乱数生成器の均等分布次元を計算する
     *
     * PIS法(原瀬)によって疑似乱数生成器の出力の均等分布次元を計算する
     * アルゴリズム
     *
     * @tparam T 疑似乱数生成器の出力の型
     */
    template<typename T> class AlgorithmEquidsitribution {

        /**
         * GF(2)ベクトルとしての疑似乱数生成器
         */
        typedef linear_generator_vector<T> linear_vec;

        /*
         * 均等分布次元計算可能な疑似乱数生成器
         */
        typedef EquidistributionCalculatable<T> ECGenerator;

    public:

        /**
         * コンストラクタ
         *
         * @param rand 均等分布次元計算可能な疑似乱数生成器
         * @param bit_len_ 疑似乱数生成器の出力のビット長
         */
        AlgorithmEquidsitribution(const ECGenerator& rand, int bit_len_) {
            bit_len = bit_len_;
            size = bit_len + 1;
            basis = new linear_vec * [size];
            stateBitSize = rand.bitSize();
            for (int i = 0; i < bit_len; i++) {
                basis[i] = new linear_vec(rand, i);
            }
            basis[bit_len] = new linear_vec(rand);
            basis[bit_len]->next_state(bit_len);
        }

        /**
         * デストラクタ
         */
        ~AlgorithmEquidsitribution() {
            for (int i = 0; i < size; i++) {
                delete basis[i];
            }
            delete[] basis;
        }

        int get_all_equidist(int veq[]);
        int get_equidist(int *sum_equidist);
    private:
        int get_equidist_main(int bit_len);
        void adjust(int new_len);

        /**
         * 標準基底+1個のベクトルからなる配列
         */
        linear_vec **basis;

        /**
         * 疑似乱数生成器の出力のビット長
         */
        int bit_len;

        /**
         * 疑似乱数生成器の状態空間のビット数
         */
        int stateBitSize;

        /**
         * basis の配列の要素数
         */
        int size;
    };

    /**
     * ビット長の調整と主要項の調整
     *
     * @param new_len 新しいビット長
     */
    template<typename T>
    void AlgorithmEquidsitribution<T>::adjust(int new_len) {
        using namespace std;

        T mask = (~static_cast<T>(0)) << (sizeof(T) * 8 - new_len);
        for (int i = 0; i < size; i++) {
            basis[i]->next = basis[i]->next & mask;
            if (basis[i]->next == 0) {
                basis[i]->next_state(new_len);
            }
        }
    }

#if defined(DEBUG)
    /**
     * デバッグ出力
     */
    template<typename T>
    void linear_generator_vector<T>::debug_print() {
        using namespace std;

        cout << "debug ====" << endl;
        cout << "count = " << dec << count << endl;
        cout << "zero = " << zero << endl;
        cout << "next = " << hex << next << endl;
        cout << "debug ====" << endl;
        //rand->debug_print();
    }
#else
    template<typename T>
    void linear_generator_vector<T>::debug_print() {
    }
#endif

    /**
     * vビット精度の均等分布次元を計算する。
     * v = 1 から \b bit_len までの均等分布次元を計算して、\b veq[]
     * に入れる。返却値はv=1からbit_len までの均等分布次元の理論的上限との
     * 差の総和である。
     *
     * \b 注意： AlgorithmEquidistribution をコンストラクトしてから、
     * get_all_dquidist() または、get_equidist() のどちらか一方しか
     * 呼び出すことはできない。
     *
     * @param[out] veq v ビット精度の均等分布次元の配列
     * @return 実際のvビット精度の均等分布次元と理論的上限の差の総和
     */
    template<typename T>
    int AlgorithmEquidsitribution<T>::get_all_equidist(int veq[]) {
        using namespace std;

        int sum = 0;

        veq[bit_len - 1] = get_equidist_main(bit_len);
#if defined(DEBUG)
        for (int i = 0; i < size; i++) {
            basis[i]->debug_print();
        }
#endif
        sum += stateBitSize / bit_len - veq[bit_len - 1];
        bit_len--;
        for (; bit_len >= 1; bit_len--) {
            adjust(bit_len);
            veq[bit_len - 1] = get_equidist_main(bit_len);
            sum += stateBitSize / bit_len - veq[bit_len - 1];
        }
        return sum;
    }

    /**
     * \b bit_len ビット精度の均等分布次元を計算する。
     * sum_equidist には、1 から bit_len -1 までの均等分布次元と理論的上限の
     * 差の総和が返される。
     *
     * @param sum_equidist 1からbit_len -1 までの理論的上限との差の総和
     * @return \b bit_len ビット精度の均等分布次元
     */
    template<typename T>
    int AlgorithmEquidsitribution<T>::get_equidist(int *sum_equidist) {
        using namespace std;

        int veq = get_equidist_main(bit_len);
        int sum = 0;
        bit_len--;
        for (; bit_len >= 1; bit_len--) {
            adjust(bit_len);
            sum += stateBitSize / bit_len - get_equidist_main(bit_len);
        }
        *sum_equidist = sum;
        return veq;
    }

    /**
     * ベクトルの加法
     * @param src このベクトルに足す相手のベクトル
     */
    template<typename T>
    void linear_generator_vector<T>::add(
        const linear_generator_vector<T>& src) {
        using namespace std;

        rand->add(*src.rand);
        next ^= src.next;
    }

    /**
     * 疑似乱数生成器の状態遷移
     *
     * @param bit_len MSB からの bit 長
     */
    template<typename T>
    void linear_generator_vector<T>::next_state(int bit_len) {
        using namespace std;

        if (zero) {
            return;
        }
        int zero_count = 0;
        next = rand->generate(bit_len);
        count++;
        while (next == 0) {
            zero_count++;
            if (zero_count > rand->bitSize() * 2) {
                zero = true;
                if (rand->isZero()) {
                    zero = true;
                }
                break;
            }
            next = rand->generate(bit_len);
            count++;
        }
    }

    /**
     * PIS法によるvビット精度均等分布次元の計算
     *
     * @param v MSBからのビット長
     * @return v ビット精度均等分布次元
     */
    template<typename T>
    int AlgorithmEquidsitribution<T>::get_equidist_main(int v) {
        using namespace std;
        using namespace NTL;
        int bit_len = v;
        int pivot_index;
        int old_pivot = 0;

        pivot_index = calc_1pos(basis[bit_len]->next);
        while (!basis[bit_len]->zero) {
#if defined(DEBUG)
            if (pivot_index != calc_1pos(basis[pivot_index]->next)) {
                cerr << "pivot error 1" << endl;
                cerr << "pivot_index:" << dec << pivot_index << endl;
                cerr << "calc_1pos:" << dec
                     << calc_1pos(basis[pivot_index]->next) << endl;
                cerr << "next:" << hex << basis[pivot_index]->next << endl;
                throw new std::logic_error("pivot error 1");
            }
#endif
            // アルゴリズムとして、全部のcount を平均的に大きくしたい。
            // 従って count の小さい方を変化させたい
            if (basis[bit_len]->count > basis[pivot_index]->count) {
                swap(basis[bit_len], basis[pivot_index]);
            }
            basis[bit_len]->add(*basis[pivot_index]);
            // add の結果 next の最後の1 は必ず 0 になる。
            // 全部0なら次の状態に進める。（内部でcount が大きくなる）
            if (basis[bit_len]->next == 0) {
                basis[bit_len]->next_state(bit_len);
                pivot_index = calc_1pos(basis[bit_len]->next);

            // 全部0でなければ、pivot_index は小さくなる。
            // pivot_index が 0 になれば最上位bit のみ1なので
            // 次の add で全部0になる。
            } else {
                old_pivot = pivot_index;
                pivot_index = calc_1pos(basis[bit_len]->next);
                if (old_pivot <= pivot_index) {
                    cerr << "pivot error 2" << endl;
                    throw new std::logic_error("pivot error 2");
                }
            }
        }

        // 計算終了したので最長のベクトルを求める。（長いとはcountが少ないこと）
        int min_count = basis[0]->count;
        for (int i = 1; i < bit_len; i++) {
            if (min_count > basis[i]->count) {
                min_count = basis[i]->count;
            }
        }
        if (min_count > stateBitSize / bit_len) {
            cout << "over theoretical bound " << bit_len << endl;
            for(int i = 0; i < size; i++) {
                basis[i]->debug_print();
            }
            throw new std::logic_error("over theoretical bound");
        }
        return min_count;
    }
}
#endif // MTTOOLBOX_ALGORITHM_EQUIDISTRIBUTION_HPP

