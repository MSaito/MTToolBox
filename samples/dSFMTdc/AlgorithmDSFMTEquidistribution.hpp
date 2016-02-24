#pragma once
#ifndef MTTOOLBOX_ALGORITHM_DSFMT_EQUIDISTRIBUTION_HPP
#define MTTOOLBOX_ALGORITHM_DSFMT_EQUIDISTRIBUTION_HPP
/**
 * @file AlgorithmDSFMTEquidistribution.hpp
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
 * Copyright (C) 2015 Mutsuo Saito, Makoto Matsumoto, Manieth Corp
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
#include <MTToolBox/util.hpp>
#include <limits.h>

namespace MTToolBox {
    struct DSFMTInfo {
        int bitSize;
        int elementNo;
    };

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
    template<typename U, typename SIMDGenerator>
    class dsfmt_linear_generator_vector {
    public:

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
        dsfmt_linear_generator_vector<U, SIMDGenerator>(
            const SIMDGenerator& generator,
            DSFMTInfo& info)
            {
#if defined(MTTOOLBOX_USE_TR1)
            std::tr1::shared_ptr<SIMDGenerator> r(new SIMDGenerator(generator));
#else
            std::shared_ptr<SIMDGenerator> r(new SIMDGenerator(generator));
#endif
            rand = r;
            count = 0;
            zero = false;
            setZero(next);
            this->info = info;
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
        dsfmt_linear_generator_vector<U, SIMDGenerator>(
            const SIMDGenerator& generator,
            int bit_pos, DSFMTInfo& info) {
#if defined(MTTOOLBOX_USE_TR1)
            std::tr1::shared_ptr<SIMDGenerator> r(new SIMDGenerator(generator));
#else
            std::shared_ptr<SIMDGenerator> r(new SIMDGenerator(generator));
#endif
            rand = r;
            rand->setZero();
            count = 0;
            zero = false;
            next = getOne<U>() << (bit_size<U>() - bit_pos - 1);
            this->info = info;
#if defined(DEBUG)
            cout << "DEBUG:" << dec << bit_pos << ":"
                 << hex << next << endl;
#endif
        }

        void add(const dsfmt_linear_generator_vector<U, SIMDGenerator>& src);
        void get_next(int bit_len);
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
#if defined(MTTOOLBOX_USE_TR1)
        std::tr1::shared_ptr<SIMDGenerator> rand;
#else
        std::shared_ptr<SIMDGenerator> rand;
#endif
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
        DSFMTInfo info;
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
     *\endenglish
     */
    template<typename U, typename SIMDGenerator>
    class AlgorithmDSFMTEquidistribution {

        /**
         *\japanese
         * GF(2)ベクトルとしての疑似乱数生成器
         *\endjapanese
         *
         *\english
         * Pseudo random number generator as a vector.
         *\endenglish
         */
        typedef dsfmt_linear_generator_vector<U, SIMDGenerator> linear_vec;
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
         *\endenglish
         */
        AlgorithmDSFMTEquidistribution(const SIMDGenerator& rand,
                                      int bit_length,
                                      DSFMTInfo info,
                                      int maxbitsize) {
#if defined(DEBUG)
            cout << "AlgorithmDSFMTEquidistribution constructer start" << endl;
            cout << "bit_length = " << dec << bit_length << endl;
#endif
            int bit_len = bit_length;
            int bit_size = bit_len * info.elementNo;
#if defined(DEBUG)
            cout << "bit_size = " << dec << bit_size << endl;
#endif
            this->info = info;
            size = bit_size + 1;
            basis = new linear_vec * [size];
            stateBitSize = maxbitsize;
            for (int i = 0; i < bit_size; i++) {
                basis[i] = new linear_vec(rand, i, info);
            }
            basis[bit_size] = new linear_vec(rand, info);
            basis[bit_size]->next_state(bit_len);
#if defined(DEBUG)
            cout << "zero = " << dec << basis[bit_size]->zero << endl;
            cout << "count = " << dec << basis[bit_size]->count << endl;
            cout << "AlgorithmDSFMTEquidistribution constructer end" << endl;
#endif
        }

        /**
         *\japanese
         * デストラクタ
         *\endjapanese
         *\english
         * Destructor
         *\endenglish
         */
        ~AlgorithmDSFMTEquidistribution() {
            for (int i = 0; i < size; i++) {
                delete basis[i];
            }
            delete[] basis;
        }

        int get_equidist(int bitLen);
    private:
        int get_equidist_main(int bit_len);

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
        //int bit_len;

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
        DSFMTInfo info;
        int size;
    };

#if defined(DEBUG)
    /**
     *\japanese
     * デバッグ出力
     *\endjapanese
     *\english
     * debug output
     *\endenglish
     */
    template<typename U, typename V>
    void dsfmt_linear_generator_vector<U, V>::debug_print() {
        using namespace std;

        cout << "count = " << dec << count;
        cout << " zero = " << dec << zero;
        cout << " next = " << hex << next << endl;
    }
#else
    template<typename U, typename V>
    void dsfmt_linear_generator_vector<U, V>::debug_print() {
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
    template<typename U, typename V>
    int AlgorithmDSFMTEquidistribution<U, V>::get_equidist(int bitLen)
    {
        using namespace std;
        return get_equidist_main(bitLen);
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
    template<typename U, typename SIMDGenerator>
    void dsfmt_linear_generator_vector<U, SIMDGenerator>::add(
        const dsfmt_linear_generator_vector<U, SIMDGenerator>& src) {
        using namespace std;

        rand->add(*src.rand);
        next ^= src.next;
    }

    /**
     * nextを作る。
     *
     */
    template<typename U, typename SIMDGenerator>
    void dsfmt_linear_generator_vector<U, SIMDGenerator>::get_next(int bitSize) {
        using namespace std;
        U w = rand->generate();
#if defined(DEBUG) && 0
        cout << "w = " << w << endl;
#endif
        bitSize = bitSize * info.elementNo;
        setZero(next);
        int k = info.bitSize - 1;
        for (int i = 0; i < info.elementNo; i++) {
            uint64_t mask = UINT64_C(0x0008000000000000);
            for (int j = 0; j < bitSize; j += info.elementNo) {
                if (w.u64[i] & mask) {
                    setBitOfPos(&next, k, 1);
                } else {
                    setBitOfPos(&next, k, 0);
                }
                k--;
                mask = mask >> 1;
            }
        }
#if defined(DEBUG)
        if (!isZero(next)) {
            cout << "bitSize = " << dec << bitSize;
            cout << " info.elementNo = " << dec << info.elementNo << endl;
            cout << "get_next w = " << hex << w << endl;
            cout << "get_next next = " << next << endl;
        }
#endif
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
    template<typename U, typename V>
    void dsfmt_linear_generator_vector<U, V>::next_state(int bitLen) {
        using namespace std;

        if (zero) {
            return;
        }
        int zero_count = 0;
        get_next(bitLen);
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
            get_next(bitLen);
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
    template<typename U, typename SIMDGenerator>
    int AlgorithmDSFMTEquidistribution<U, SIMDGenerator>::
    get_equidist_main(int v)
    {
        using namespace std;
        using namespace NTL;
#if defined(DEBUG)
        cout << "get_equidist_main start" << endl;
        cout << "v = " << dec << v << endl;
#endif
        int bitSize = v * info.elementNo;
        int pivot_index;
        int old_pivot = 0;

        pivot_index = calc_1pos(basis[bitSize]->next);
#if defined(DEBUG)
        cout << "get_equidist_main step 1" << endl;
#endif
        while (!basis[bitSize]->zero) {
#if defined(DEBUG)
            cout << "get_equidist_main step 2" << endl;
#endif
#if defined(DEBUG)
            if (pivot_index == -1) {
                cout << "pivot_index = " << dec << pivot_index << endl;
                cout << "zero = " << basis[bitSize]->zero << endl;
                cout << "next = " << hex << basis[bitSize]->next << endl;
                throw new std::logic_error("pivot error 0");
            }
            if (pivot_index >= bitSize) {
                cout << "pivot_index = " << dec << pivot_index << endl;
                cout << "bitSize = " << bitSize << endl;
                cout << "next = " << hex << basis[bitSize]->next << endl;
                throw new std::logic_error("pivot error 0.1");
            }
            if (pivot_index != calc_1pos(basis[pivot_index]->next)) {
                cerr << "pivot error 1" << endl;
                cerr << "pivot_index:" << dec << pivot_index << endl;
                cerr << "calc_1pos:" << dec
                     << calc_1pos(basis[pivot_index]->next) << endl;
                cerr << "next:" << hex << basis[pivot_index]->next << endl;
                for (int i = 0; i < bitSize; i++) {
                    cerr << dec << i << ":" << hex << basis[i]->next << endl;
                }
                throw new std::logic_error("pivot error 1");
            }
#endif
            // アルゴリズムとして、全部のcount を平均的に大きくしたい。
            // 従って count の小さい方を変化させたい
            if (basis[bitSize]->count > basis[pivot_index]->count) {
                swap(basis[bitSize], basis[pivot_index]);
            }
#if defined(DEBUG)
            cout << "before add bitSize next = " << hex
                 << basis[bitSize]->next << endl;
            cout << "before add pivot   next = " << hex
                 << basis[pivot_index]->next << endl;
#endif
            basis[bitSize]->add(*basis[pivot_index]);
#if defined(DEBUG)
            cout << "after  add bitSize next = " << hex
                 << basis[bitSize]->next << endl;
#endif
            // add の結果 next の最後の1 は必ず 0 になる。
            // 全部0なら次の状態に進める。（内部でcount が大きくなる）
            if (isZero(basis[bitSize]->next)) {
                basis[bitSize]->next_state(v);
                pivot_index = calc_1pos(basis[bitSize]->next);
#if defined(DEBUG)
                cout << "zero" << endl;
                cout << "pivot_index = " << dec << pivot_index << endl;
                if (pivot_index >= bitSize) {
                    cout << "pivot_index = " << dec << pivot_index << endl;
                    cout << "bitSize = " << bitSize << endl;
                    cout << "next = " << hex << basis[bitSize]->next << endl;
                    throw new std::logic_error("pivot error 1.1");
                }
                if (basis[bitSize]->zero) {
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
                pivot_index = calc_1pos(basis[bitSize]->next);
                if (pivot_index >= bitSize) {
                    cout << "pivot_index = " << dec << pivot_index << endl;
                    cout << "bitSize = " << bitSize << endl;
                    cout << "next = " << hex << basis[bitSize]->next << endl;
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
#if defined(DEBUG)
        for (int i = 0; i < bitSize; i++) {
            cout << dec << i << ": count = " << basis[i]->count << endl;
        }
#endif
        int min_count = basis[0]->count;
//        int min_count = INT_MAX;
        for (int i = 1; i < bitSize; i++) {
            if (basis[i]->zero) {
                continue;
            }
            if (min_count > basis[i]->count) {
                min_count = basis[i]->count;
            }
        }
        int result = min_count;
        if (result > stateBitSize / bitSize) {
//        if (result > stateBitSize / (bitSize / info.elementNo)) {
            cerr << "over theoretical bound" << endl;
            cout << "bitSize = " << dec << bitSize << endl;
            cout << "min_count = " << dec << min_count << endl;
            cout << "stateBitSize = " << dec << stateBitSize << endl;
            cout << "elementNo = " << dec << info.elementNo << endl;
            cout << "result = " << dec << result << endl;
            cerr << basis[0]->rand->getParamString() << endl;
#if 0
            for(int i = 0; i < size; i++) {
                basis[i]->debug_print();
            }
#endif
            throw new std::logic_error("over theoretical bound");
        }
#if defined(DEBUG)
        cout << "get_equidist_main end" << endl;
#endif
        return result;
    }

    template<typename U, typename SIMDGenerator>
    int calc_dSFMT_equidistribution(const SIMDGenerator& rand,
                                   int veq[],
                                   int bit_len,
                                   DSFMTInfo& info,
                                   int mexp)
    {
        int state_inc = 1;
        int weight_max = info.elementNo;
        int state_max = weight_max;
        int weight_dec = state_inc;
        for (int i = 0; i < bit_len; i++) {
            veq[i] = INT_MAX;
        }
        int veq_weight[bit_len];
        for (int sm = 0; sm < state_max; sm += state_inc) {
            for (int i = 0; i < bit_len; i++) {
                veq_weight[i] = -1;
            }
            for (int wm = weight_dec; wm <= weight_max; wm += weight_dec) {
#if 0
                cout << "start_mode = " << dec << sm;
                cout << " weight_mode = " << dec << wm << endl;
#endif
                //SIMDGenerator work = rand;
                //work.setStartMode(sm);
                //work.setWeightMode(wm);
                // previous set
                //work.generate();
                for (int v = 1; v <= bit_len; v++) {
                    SIMDGenerator work = rand;
                    work.setStartMode(sm);
                    work.setWeightMode(wm);
                    // previous set
                    work.generate();
                    AlgorithmDSFMTEquidistribution<U, SIMDGenerator>
                        ase(work, v, info, rand.bitSize());
                    int e = ase.get_equidist(v);
#if 0
                    cout << "min_count = " << dec << e;
#endif
                    e = e * info.elementNo - (info.elementNo - wm);
                    if (e > mexp / v) {
                        cerr << "over theoretical bound" << endl;
                        cout << "start_mode = " << dec << sm;
                        cout << " weight_mode = " << dec << wm;
                        cout << " mexp = " << dec << mexp;
                        cout << " e = " << dec << e;
                        cout << " v = " << dec << v << endl;
                        throw new std::logic_error("over theoretical bound");
                    }
#if 0
                    cout << "\tk(" << dec << v << ") = " << dec
                         << e << endl;
#endif
                    // max
                    if (e > veq_weight[v - 1]) {
                        veq_weight[v - 1] = e;
                    }
                }
            }
            for (int i = 0; i < bit_len; i++) {
                // min
                if (veq[i] > veq_weight[i]) {
                    veq[i] = veq_weight[i];
                }
            }
        }
        int sum = 0;
        for (int i = 1; i <= bit_len; i++) {
            sum += mexp / i - veq[i - 1];
        }
        return sum;
    }

}
#if defined(MTTOOLBOX_USE_TR1)
#undef MTTOOLBOX_USE_TR1
#endif
#endif // MTTOOLBOX_ALGORITHM_EQUIDISTRIBUTION_HPP
