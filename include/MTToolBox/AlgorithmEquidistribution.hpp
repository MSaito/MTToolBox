#ifndef MTTOOLBOX_ALGORITHM_EQUIDISTRIBUTION_HPP
#define MTTOOLBOX_ALGORITHM_EQUIDISTRIBUTION_HPP
/**
 * @file AlgorithmEquidistribution.hpp
 *
 *\japanese
 * @brief 疑似乱数生成器の出力の均等分布次元を計算する。
 *
 * PIS法[1](原瀬)によって疑似乱数生成器の出力の均等分布次元を計算するアルゴリズム
 *\endjapanese
 *
 *\english
 * @brief Calculate dimension of equi-distribution of output of pseudo
 * random number generators.
 *
 * Algorithm that calculates dimension of equi-distribution of output
 * of pseudo random number generators using PIS method[1](S. Harase).
 *\endenglish
 *
 * @author Mutsuo Saito (Manieth Corp.)
 * @author Makoto Matsumoto (Hiroshima University)
 *
 * Copyright (C) 2013, 2016 Mutsuo Saito, Makoto Matsumoto
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
#if __cplusplus >= 201103L
#include <memory>
#else
#if defined(__clang__)
#if __has_include(<tr1/memory>)
#define MTTOOLBOX_USE_TR1
#include <tr1/memory>
#else
#include <memory>
#endif // __has_indlude
#else  // not clang
#define MTTOOLBOX_USE_TR1
#include <tr1/memory>
#endif // clang
#endif // cplusplus version
#include <stdexcept>
#include <MTToolBox/EquidistributionCalculatable.hpp>
#include <MTToolBox/period.hpp>
#include <MTToolBox/util.hpp>

namespace MTToolBox {
#if defined(MTTOOLBOX_USE_TR1)
    using std::tr1::shared_ptr;
#else
    using std::shared_ptr;
#endif
    /**
     * @class linear_generator_vector
     *\japanese
     * @brief GF(2)ベクトルとしてのGF(2)疑似乱数生成器
     *
     * このクラスはGF(2)ベクトルとしての疑似乱数生成器を表す。また
     * GF(2)係数多項式をv個まとめたベクトルとしても扱うが、多項式特有の
     * 処理としては不定元倍(x をかける）ことのみをサポートする。これは
     * next_state()メソッドによって行われる。
     *
     * @tparam U 疑似乱数生成器の出力の型, 符号なし型でなければならない
     *\endjapanese
     *
     *\english
     * @brief GF(2) pseudo random number generator as a GF(2) vector.
     *
     * This class treats a GF(2) pseudo random number generator as a
     * GF(2) vector. And also the class treats the generator as a
     * vector of polynomials over GF(2).
     *
     * @tparam U type of output of pseudo random number
     * generator. Should be unsigned type.
     *\endenglish
     */
    template<typename U>
    class linear_generator_vector {
    public:

        /**
         *\japanese
         * 均等分布次元計算可能な疑似乱数生成器
         *\endjapanese
         *
         *\english
         * pseudo random number generator which can calculate
         * dimension of equi-distribution.
         *\endenglish
         */
        typedef EquidistributionCalculatable<U> ECGenerator;

        /**
         *\japanese
         * コンストラクタ
         *
         * @param generator 均等分布次元を計算するGF(2)疑似乱数生成器
         *\endjapanese
         *
         * Constructor
         * @param generator GF(2) pseudo random number generator,
         * whose dimension of equi-distribution will be calculated.
         *\english
         *\endenglish
         */
        linear_generator_vector<U>(const ECGenerator& generator) {
            shared_ptr<ECGenerator> r(generator.clone());
            rand = r;
            //rand->seed(1);
            count = 0;
            zero = false;
            setZero(next);
        }

        /**
         *\japanese
         * 標準基底コンストラクタ
         *
         * 標準基底を構成するベクトルのひとつを生成する。
         * 出力の特定のビットに1度だけ1を出力して、その後はずっと0を出力する
         * 疑似乱数生成器のコンストラクタ。例えば、1, 0, 0, 0, .... と出力する、
         * あるいは 8, 0, 0, 0, .... と出力するなど。
         *
         * @param generator GF(2)疑似乱数生成器
         * @param bit_pos 1 の位置
         *\endjapanese
         *
         *\english
         * Constructor for standard basis.
         *
         * Construct a vector which consists standard basis.
         * As a generator, this generates a number include only one bit once,
         * and then generate zero forever. For instance, it generates
         * 1, 0, 0, 0, ... or 8, 0, 0, 0, ...
         *\endenglish
         */
        linear_generator_vector<U>(const ECGenerator& generator,
                                   int bit_pos) {
            shared_ptr<ECGenerator> r(generator.clone());
            rand = r;
            rand->setZero();
            count = 0;
            zero = false;
            next = getOne<U>() << (bit_size<U>() - bit_pos - 1);
#if defined(DEBUG)
            cout << "DEBUG:" << dec << bit_pos << ":"
                 << hex << next << endl;
#endif
        }

        void add(const linear_generator_vector<U>& src);
        void next_state(int bit_len);
        void debug_print();

        /**
         *\japanese
         * GF(2)線形疑似乱数生成器
         *\endjapanese
         *
         *\english
         * A GF(2) linear pseudo random number generator
         *\endenglish
         */
        shared_ptr<ECGenerator> rand;
        /**
         *\japanese
         * next_state() が呼ばれた回数
         * これは多項式としてみた場合の次数に関係がある。
         *\endjapanese
         *
         *\english
         * number which shows counts of next_state() called.
         * As a polynomial, this is related to degree of polynomial.
         *\endenglish
         */
        int count;

        /**
         *\japanese
         * ゼロベクトルであるかどうかを示す。
         *\endjapanese
         *
         *\english
         * Shows if zero vector or not.
         *\endenglish
         */
        bool zero;

        /**
         *\japanese
         * 疑似乱数生成器の最新の出力（の上位vビット）または
         * 多項式の最高次の係数のなすベクトル
         *\endjapanese
         *
         *\english
         * latest output of pseudo random number generator, or, a
         * GF(2) vector consists of coefficients of highest degree
         * term of polynomial.
         *\endenglish
         */
        U next;
    };

    /**
     * @class AlgorithmEquidistribution
     *\japanese
     * @brief 疑似乱数生成器の均等分布次元を計算する
     *
     * PIS法(原瀬)によって疑似乱数生成器の出力の均等分布次元を計算する
     * アルゴリズム
     *
     * @tparam U 疑似乱数生成器の出力の型
     *\endjapanese
     *
     *\english
     * @brief Calculate dimension of equi-distribution of output of
     * pseudo random number generators.
     *
     * Algorithm that calculates dimension of equi-distribution of
     * output of pseudo random number generators using PIS
     * method[1](S. Harase).
     * @tparam type of output of pseudo random number generator.
     * @tparam type of output of parameter generator.
     *\endenglish
     */
    template<typename U>
    class AlgorithmEquidistribution {

        /**
         *\japanese
         * GF(2)ベクトルとしての疑似乱数生成器
         *\endjapanese
         *
         *\english
         * Pseudo random number generator as a vector.
         *\endenglish
         */
        typedef linear_generator_vector<U> linear_vec;

        /*
         *\japanese
         * 均等分布次元計算可能な疑似乱数生成器
         *\endjapanese
         *\english
         * Pseudo random number generator which can calculate dimension
         * of equi-distribution.
         *\endenglish
         */
        typedef EquidistributionCalculatable<U> ECGenerator;

    public:

        /**
         *\japanese
         * コンストラクタ
         *
         * PIS法の特徴としてvビット精度均等分布次元を計算する際に、k(v+1)
         * の計算時の中間結果を利用してk(v)の計算の手間を省くことができる。この
         * クラスではその特徴を反映してk(v)のvをbit_len から1まで変化さ
         * せて一度に求めることができるようになっている。
         *
         * @param rand 均等分布次元計算可能な疑似乱数生成器
         * @param bit_length 均等分布次元を計算するMSBからのビット長, k(v)のv
         * の最初の値
         * @param mexp メルセンヌ指数、指定しなければ状態空間サイズ。
         * 可約生成器の場合は指定する必要がある。
         *\endjapanese
         *
         *\english
         * Constructor
         *
         * PIS method can calculates dimension of equi-distribution of
         * v-bit accuracy k(v) using intermediate result of
         * calculation of k(v+1). This class calculate multiple k(v)
         * varying \b v from \b bit_length to 1.
         *
         * @param rand pseudo random number generator
         * @param bit_length bit length from MSB to calculate dimension
         * of equi-distribution. This is first v of k(v).
         * @param mexp Mersenne Exponent. If not specified, internal state
         * size is used. If generator is reducible, this parameter should
         * be specified.
         *\endenglish
         */
        AlgorithmEquidistribution(const ECGenerator& rand, int bit_length,
                                  int mexp = 0) {
            bit_len = bit_length;
            size = bit_len + 1;
            basis = new linear_vec * [static_cast<unsigned long>(size)];
            if (mexp == 0) {
                stateBitSize = rand.bitSize();
            } else {
                stateBitSize = mexp;
            }
            for (int i = 0; i < bit_len; i++) {
                basis[i] = new linear_vec(rand, i);
            }
            basis[bit_len] = new linear_vec(rand);
            basis[bit_len]->next_state(bit_len);
        }

        /**
         *\japanese
         * デストラクタ
         *\endjapanese
         *\english
         * Destructor
         *\endenglish
         */
        ~AlgorithmEquidistribution() {
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
         *\japanese
         * 標準基底+1個のベクトルからなる配列。
         * basis という名前だが基底ではなく格子の生成集合という方が正しい。
         * k(v)を算出した時点では、ゼロベクトルとなっている
         * ベクトルを除けばv次元の格子の基底となっている。
         *\endjapanese
         *
         *\english
         * An array consists of standard basis plus one vector.
         * This array has a name \b basis, but this array is not
         * basis, this array should be called generating set of
         * lattice. But this array becomes basis of v dimensional
         * lattice when k(v) has calculated excluding zero vector.
         *\endenglish
         */
        linear_vec **basis;

        /**
         *\japanese
         * vビット精度均等分布次元の計算の v
         *\endjapanese
         *\english
         * \b v of dimension of equi-distribution with v-bit accuracy.
         *\endenglish
         */
        int bit_len;

        /**
         *\japanese
         * 疑似乱数生成器の状態空間のビット数
         *\endjapanese
         *\english
         * Number of bits in internal state of pseudo random number
         * generator.
         *\endenglish
         */
        int stateBitSize;

        /**
         *\japanese
         * basis の配列の要素数
         *\endjapanese
         *
         *\english
         * Size of array \b basis.
         *\endenglish
         */
        int size;
    };

    /**
     *\japanese
     * ビット長の調整と主要項の調整
     *
     * @param new_len 新しいビット長
     *\endjapanese
     *
     *\english
     * Adjust bit length and coefficients of primary terms.
     * This is needed to use previous intermediate result for
     * next calculation.
     *\endenglish
     */
    template<typename U>
    void AlgorithmEquidistribution<U>::adjust(int new_len) {
        using namespace std;
        U tmp;
        setZero(tmp);
        U mask = ~tmp << (bit_size<U>() - new_len);
        for (int i = 0; i < size; i++) {
            basis[i]->next = basis[i]->next & mask;
            if (isZero(basis[i]->next)) {
                basis[i]->next_state(new_len);
            }
        }
    }

#if defined(DEBUG)
    /**
     *\japanese
     * デバッグ出力
     *\endjapanese
     *\english
     * debug output
     *\endenglish
     */
    template<typename U>
    void linear_generator_vector<U>::debug_print() {
        using namespace std;

        cout << "debug ====" << endl;
        cout << "count = " << dec << count << endl;
        cout << "zero = " << zero << endl;
        cout << "next = " << hex << next << endl;
        cout << "debug ====" << endl;
    }
#else
    template<typename U>
    void linear_generator_vector<U>::debug_print() {
    }
#endif

    /**
     *\japanese
     * vビット精度の均等分布次元を計算する。
     * v = \b bit_len から 1までの均等分布次元を計算して、\b veq[]
     * に入れる。返却値はv=1からbit_len までの均等分布次元の理論的上限との
     * 差の総和である。
     *
     * \warning AlgorithmEquidistribution をコンストラクトしてから、
     * get_all_equidist() または、get_equidist() のどちらか一方を１回し
     * か呼び出すことはできない。
     *
     * @param[out] veq v ビット精度の均等分布次元の配列
     * @return 実際のvビット精度の均等分布次元と理論的上限の差の総和
     *\endjapanese
     *
     *\english
     *
     * Calculate dimension of equi-distribution with v-bit accuracy.
     *
     * Calculate dimension of equi-distribution with v-bit accuracy
     * k(v) for v = \b bit_length to 1, and set them into an array \b
     * veq[].  The return value is sum of d(v)s, which are difference
     * between k(v) and theoretical upper bound at \b v.
     *
     * \warning Only one of get_all_equidist() or get_equidist() can
     * be called only one time. (Bad interface)
     *
     * @param[out] veq an array of k(v)
     * @return sum of d(v)s
     *
     *\endenglish
     */
    template<typename U>
    int AlgorithmEquidistribution<U>::get_all_equidist(int veq[]) {
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
     *\japanese
     * vビット精度の均等分布次元を計算する。
     *
     * コンストラクタで指定したbit_length についてk(v)を計算して返す。
     * sum_equidist には、1 から bit_len -1 までの均等分布次元と理論的上限の
     * 差の総和が返される。
     *
     * \warning AlgorithmEquidistribution をコンストラクトしてから、
     * get_all_equidist() または、get_equidist() のどちらか一方を１回し
     * か呼び出すことはできない。
     *
     * @param sum_equidist 1からbit_len -1 までの理論的上限との差の総和
     * @return \b bit_length ビット精度の均等分布次元
     *\endjapanese
     *
     *\english
     * Calculate dimension of equi-distribution with v-bit accuracy.
     *
     * Calculate dimension of equi-distribution with v-bit accuracy
     * k(v) for v = \b bit_length.
     * \b sum_equidist is sum of d(v)s, which are difference
     * between k(v) and theoretical upper bound at \b v.
     *
     * \warning Only one of get_all_equidist() or get_equidist() can
     * be called only one time. (Bad interface)
     *
     * @param[out] sum_equidist sum of d(v)s
     * @return k(bit_length)
     *
     *\endenglish
     */
    template<typename U>
    int AlgorithmEquidistribution<U>::get_equidist(int *sum_equidist) {
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
     *\japanese
     * ベクトルの加法
     * @param src このベクトルに足す相手のベクトル
     *\endjapanese
     *
     *\english
     * Vector addition
     * @param src source vector to be added to this vector
     *\endenglish
     */
    template<typename U>
    void linear_generator_vector<U>::add(
        const linear_generator_vector<U>& src) {
        using namespace std;

        rand->add(*src.rand);
        next ^= src.next;
    }

    /**
     *\japanese
     * 疑似乱数生成器の状態遷移
     *
     * 多項式ベクトルとしてみると、すべての多項式を不定元倍する。
     *
     * @param bit_len MSB からの bit 長
     *\endjapanese
     *
     *\english
     * State transition of pseudo random number generator.
     *
     * As a vector of polynomial, multiply by an indeterminate.
     * @param bit_len bit length from MSB, or \b v of k(v) currently
     * calculating.
     *\endenglish
     */
    template<typename U>
    void linear_generator_vector<U>::next_state(int bit_len) {
        using namespace std;

        if (zero) {
            return;
        }
        int zero_count = 0;
        next = rand->generate(bit_len);
        count++;
        while (isZero(next)) {
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
     *\japanese
     * PIS法によるvビット精度均等分布次元の計算のメインとなるメソッド
     *
     * @param v MSBからのビット長
     * @return v ビット精度均等分布次元
     *\endjapanese
     *
     *\english
     * Main method of calculation of dimension of equi-distribution
     * with v-bit accuracy.
     *
     * @param v v of k(v)
     * @return k(v)
     *\endenglish
     */
    template<typename U>
    int AlgorithmEquidistribution<U>::get_equidist_main(int v) {
        using namespace std;
        using namespace NTL;
        int bit_len = v;
        int pivot_index;
        int old_pivot = 0;

        pivot_index = calc_1pos(basis[bit_len]->next);
        while (!basis[bit_len]->zero) {
#if defined(DEBUG)
            if (pivot_index == -1) {
                cout << "pivot_index = " << dec << pivot_index << endl;
                cout << "zero = " << basis[bit_len]->zero << endl;
                cout << "next = " << hex << basis[bit_len]->next << endl;
                throw new std::logic_error("pivot error 0");
            }
            if (pivot_index >= bit_len) {
                cout << "pivot_index = " << dec << pivot_index << endl;
                cout << "bit_len = " << bit_len << endl;
                cout << "next = " << hex << basis[bit_len]->next << endl;
                throw new std::logic_error("pivot error 0.1");
            }
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
#if defined(DEBUG)
            cout << "before add bit_len next = " << hex
                 << basis[bit_len]->next << "  count = " << dec
                 << basis[bit_len]->count << endl;
            cout << "before add pivot   next = " << hex
                 << basis[pivot_index]->next << "  count = " << dec
                 << basis[pivot_index]->count << endl;
#endif
            basis[bit_len]->add(*basis[pivot_index]);
#if defined(DEBUG)
            cout << "after  add bit_len next = " << hex
                 << basis[bit_len]->next << "  count = " << dec
                 << basis[bit_len]->count << endl;
#endif
            // add の結果 next の最後の1 は必ず 0 になる。
            // 全部0なら次の状態に進める。（内部でcount が大きくなる）
            if (isZero(basis[bit_len]->next)) {
                basis[bit_len]->next_state(bit_len);
                pivot_index = calc_1pos(basis[bit_len]->next);
#if defined(DEBUG)
                cout << "zero" << endl;
                cout << "pivot_index = " << dec << pivot_index << endl;
                if (pivot_index >= bit_len) {
                    cout << "pivot_index = " << dec << pivot_index << endl;
                    cout << "bit_len = " << bit_len << endl;
                    cout << "next = " << hex << basis[bit_len]->next << endl;
                    throw new std::logic_error("pivot error 1.1");
                }
                if (basis[bit_len]->zero) {
                    cout << "loop exit condition" << endl;
                } else if (pivot_index == -1) {
                    cerr << "pivot error 1.1" << endl;
                    throw new std::logic_error("pivot error 1.2");
                }
#endif
            // 全部0でなければ、pivot_index は小さくなる。
            // pivot_index が 0 になれば最上位bit のみ1なので
            // 次の add で全部0になる。
            } else {
                old_pivot = pivot_index;
                pivot_index = calc_1pos(basis[bit_len]->next);
                if (pivot_index >= bit_len) {
                    cout << "pivot_index = " << dec << pivot_index << endl;
                    cout << "bit_len = " << bit_len << endl;
                    cout << "next = " << hex << basis[bit_len]->next << endl;
                    throw new std::logic_error("pivot error 2");
                }
                if (old_pivot <= pivot_index) {
                    cerr << "pivot error 2" << endl;
                    cerr << "old_pivot = " << dec << old_pivot << endl;
                    cerr << "pivot_index = " << dec << pivot_index << endl;
                    throw new std::logic_error("pivot error 2.1");
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
#if defined(DEBUG)
        cout << "min_count = " << min_count << endl;
        cout << "stateBitSize = " << stateBitSize << endl;
        cout << "bit_len = " << bit_len << endl;
        cout << "theoretical bound " << (stateBitSize / bit_len) << endl;
#endif
        if (min_count > stateBitSize / bit_len) {
            cerr << basis[0]->rand->getParamString() << endl;
            for(int i = 0; i < size; i++) {
                basis[i]->debug_print();
            }
            cerr << "min_count = " << min_count << endl;
            cerr << "stateBitSize = " << stateBitSize << endl;
            cerr << "bit_len = " << bit_len << endl;
            cerr << "over theoretical bound " << (stateBitSize / bit_len)
                 << endl;
            throw new std::logic_error("over theoretical bound");
        }
        return min_count;
    }
}
#if defined(MTTOOLBOX_USE_TR1)
#undef MTTOOLBOX_USE_TR1
#endif
#endif // MTTOOLBOX_ALGORITHM_EQUIDISTRIBUTION_HPP
