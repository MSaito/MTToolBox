#ifndef MTTOOLBOX_ALGORITHM_CALCULATE_PARITY_HPP
#define MTTOOLBOX_ALGORITHM_CALCULATE_PARITY_HPP
/**
 * @file AlgorithmCalculateParity.hpp
 *
 *\japanese
 * @brief 可約ジェネレータのパリティチェックベクトルを求める。
 *
 * 状態遷移関数の特性多項式hを f * q という多項式の積で表すことができて、
 * かつ f の次数が状態空間が許す限り大きなメルセンヌ指数である場合、状
 * 態遷移の周期が2^(deg(f))-1のゼロでない倍数となるように保証するような
 * ベクトルを求めることができる。ケーリー・ハミルトンの式によって、hは
 * 状態空間をannihilateするが、f, q も状態空間の部分空間をannihilateす
 * る。このクラスでは、q によってannihilateされる部分空間の基底を求める
 * ことによって、周期保証ベクトル（パリティベクトル）を求める。
 *
 * 状態空間の一部をパリティチェック用の部分空間Pとする。この空間のサイ
 * ズは状態空間を表す配列の要素のサイズになるようにする。典型的には32ビッ
 * トまたは64ビットである。状態空間からパリティチェック用の空間への射影
 * をPrとする。qのカーネルの基底を具体的に計算して、求めておく(Ker_q)。
 * Pr(Ker_q)はPの部分空間の基底になるので、PにおけるPr(Ker_q)の直交補空
 * 間の基底を計算する。その基底からひとつのベクトルを取り出してパリティ
 * チェックベクトルとする。
 *
 * 初期化後に、Pのビット列を取り出してパリティチェックベクトルとの内積
 * を計算する。内積が1ならば、Pr(Ker_q)に入っていないので、Ker_qにも入っ
 * ていない。状態空間はKer_q と Ker_f の直和なので、Ker_q に入っていな
 * ければ、Ker_f のベクトルを和の一部として含む。これにより周期が保証さ
 * れる。内積が0の場合は、Ker_q に入っているかも知れないし、入っていな
 * いかも知れないが、1ビット変更して内積が1になるように変更すれば周期が
 * 保証できる。
 *
 *\endjapanese
 *
 *\english
 * @brief Calculate the parity check vector of reducible generator.
 *
 * Let h be a characteristic polynomial of state transition function
 * of given reducible generator. h is factorized h = f * q, where
 * f is large irreducible polynomial with degree of Mersenne Exponet,
 * and q is small polynomial which may be reducible.
 * h annihilates the internal state of the generator. f and q also
 * annihilate subspaces of the internal state. By calculating
 * basis of the subspace annihilated by q, we can get the
 * period certification vector (parity vector).
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
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <iostream>
#include <inttypes.h>
#include <NTL/GF2X.h>
#include <NTL/vec_GF2.h>
#include <NTL/mat_GF2.h>
#include <MTToolBox/ReducibleGenerator.hpp>
#include <MTToolBox/util.hpp>

namespace MTToolBox {
    using namespace std;

    /**
     *\japanese
     * @brief 可約ジェネレータのパリティチェックベクトルを求める。
     *
     *\endjapanese
     *
     *\english
     * @brief Calculate the parity check vector of reducible generator.
     *\endenglish
     */
    template<typename U, typename G>
    class AlgorithmCalculateParity {

    public:
        /**
         *\japanese
         * qによってannihilateされる空間の基底を求め、パリティチェック用の定数を求める
         * @param[in] g 可約ジェネレータ
         * @param[in] f mexp次数の既約因子
         * @return パリティチェック用の定数
         *\endjapanese
         *
         *\english
         * calculate the basis of subspace annihilated by q, and get the
         * period certification vector (parity vector).
         * @param[in] g reducible generator
         * @param[in] f large irreducible factor of characteristic polynomial.
         * @return the period certification vector (parity vector)
         *\endenglish
         */
        U searchParity (G& g, const NTL::GF2X& f) {
            int mexp = g.getMexp();
#if defined(DEBUG)
            cout << "searchParity start" << endl;
            cout << "degree of f = " << deg(f) << endl;
#endif
            int maxdegree = g.bitSize();
            word_width = bit_size<U>();
            int base_num = maxdegree - mexp;
            internal_state bases[word_width];
            internal_state work_base;
            int count;
            int bit_pos = 0;

            for (int i = 0; i < word_width; i++) {
                bases[i].rg = new G(g);
                bases[i].rg->setZero();
                bases[i].zero = true;
                setZero(bases[i].next);
            }
            work_base.rg = new G(g);
            work_base.rg->setZero();
            work_base.zero = true;
            setZero(work_base.next);
#if defined(DEBUG)
            cout << "searchParity before while loop" << endl;
            cout << "word_width = " << word_width << endl;
            cout << "bit_pos = " << bit_pos << endl;
            cout << "maxdegree = " << maxdegree << endl;
            cout << "base_num = " << base_num << endl;
#endif
            while(bit_pos < maxdegree) {
                calc_basis(work_base, f, &bit_pos);
                addBase(bases, word_width, work_base);
                count = 0;
                for (int i = 0; i < word_width; i++) {
                    if (!isZero(bases[i].next)) {
                        count++;
                    }
                }
#if defined(DEBUG)
                cout << "in while loop count = " << count << endl;
#endif
                if (count >= base_num) {
                    break;
                }
            }
#if defined(DEBUG)
            cout << "searchParity after while loop" << endl;
#endif
            for (int i = 0; i < word_width; i++) {
                delete bases[i].rg;
            }
#if defined(DEBUG)
            cout << "----" << endl;
            for (int i = 0; i < word_width; i++) {
                if (isZero(bases[i].next)) {
                        continue;
                }
                int w = word_width / 4;
                cout << setw(2) << dec << i << " ";
                cout << hex << setfill('0') << setw(w) << bases[i].next << endl;
                cout << dec;
            }
            cout << "----" << endl;
#endif
            U parity = search_parity_check_vector(bases, base_num);
            g.setParityValue(parity);
            return parity;
        }

    private:
        /* internal state */
        struct internal_state {
            bool zero;
            U next;
            G * rg;
        };

        int word_width;

        /**
         *\japanese
         * q によって殲滅される部分空間の基底を求める。
         * @param[in/out] st 基底を構成するひとつのベクトル
         * @param[in] f mexp次数の既約因子
         * @param[in/out] bit_pos 状態空間内のビット位置
         *\endjapanese
         *
         *\english
         *\endenglish
         */
        void calc_basis(internal_state& st, const NTL::GF2X& f, int *bit_pos) {
#if defined(DEBUG)
            cout << "calc_basis start bit_pos = " << dec << *bit_pos << endl;
#endif
            int maxdegree = st.rg->bitSize();
            for (;*bit_pos < maxdegree;) {
#if defined(DEBUG)
                cout << "in calc_basis bit_pos = " << dec << *bit_pos << endl;
#endif
                // 全状態空間の中で1ビットだけ1を立てる
                st.rg->setOneBit(*bit_pos);
#if defined(DEBUG)
                if (st.rg->isZero()) {
                    cout << "ERROR rg is ZERO" << endl;
                }
#endif
                (*bit_pos)++;
                // fによる写像でfのカーネルの像を0にする
                annihilate(st.rg, f);
#if defined(DEBUG)
                if (st.rg->isZero()) {
                    cout << "ZERO after annihilate" << endl;
                    cout << "deg(f) = " << deg(f) << endl;
                } else {
                    cout << "NON ZERO after annihilate" << endl;
                }
#endif
                set_state(st);
                if (!st.zero) {
                    break;
                }
            }
#if defined(DEBUG)
            cout << "calc_basis end bit_pos = " << dec << *bit_pos << endl;
#endif
        }

        /**
         *\japanese
         * st.next にgenerator の次の非ゼロ値をセットする
         * @param[in/out] st 内部状態
         *\endjapanese
         *
         *\english
         * set the next non-zero value of the generator to st.next.
         * @param[in/out] st representation of internal state
         *\endenglish
         */
        void set_state(internal_state& st) {
#if defined(DEBUG)
            cout << "set_state start" << endl;
#endif
            if (st.rg->isZero()) {
                st.zero = true;
                setZero(st.next);
#if defined(DEBUG)
                cout << "set_state end zero" << endl;
#endif
                return;
            }
            st.zero = false;
            st.rg->generate();
            st.next = st.rg->getParityValue();
            while (isZero(st.next)) {
                if (st.rg->isZero()) {
                    st.zero = true;
                    break;
                }
                st.rg->generate();
                st.next = st.rg->getParityValue();
            }
#if defined(DEBUG)
            cout << "set_state end parity = "
                 << hex << setfill('0') << setw(8) << st.next << dec << endl;
#endif
        }

        /**
         *\japanese
         * GF(2)ベクトルとして状態空間を加える
         * @param[in/out] dist 加える先
         * @param[in] src 加える元
         *\endjapanese
         *
         *\english
         * add internal_state as GF(2)-vectors.
         * @param[in/out] dist distination
         * @param[in] src source
         *\endenglish
         */
        void add_state(internal_state& dist, const internal_state& src) {
            dist.rg->add(*src.rg);
            dist.next ^= src.next;
        }

        /**
         *\japanese
         * 状態空間がゼロでなければ、次状態に進める
         * @param[in/out] st 状態空間
         *\endjapanese
         *
         *\english
         * transit to next state if intenal state is not zero.
         * @param[in/out] st internal state
         *\endenglish
         */
        void get_next_state(internal_state& st) {
            if (st.zero) {
                return;
            }
            set_state(st);
        }

        /**
         *\japanese
         * 基底となるべき集合に、基底候補を加える
         * Mulders and Storjohann アルゴリズムの変種
         * @param[in/out] bases 基底
         * @param[in] size bases の大きさ
         * @param[in/out] work 基底に追加される候補となるベクトル
         *\endjapanese
         *
         *\english
         * add a candidate vector to set of vectors which will be
         * basis in the end.
         * variant of Mulders and Storjohann Algorithm
         * @param[in/out] bases basis
         * @param[in] size size of bases
         * @param[in/out] work a candidate vector to be a member of basis
         *\endenglish
         */
        void addBase(internal_state bases[], int size, internal_state& work) {
#if defined(DEBUG)
            cout << "addBase start" << endl;
#endif
            int count = bit_size<U>() * 10;
            for (;count >= 0;) {
                count--;
#if defined(DEBUG)
                cout << "addBase work.next = " << hex << work.next << endl;
#endif
                if (isZero(work.next)) {
                    get_next_state(work);
                    if (work.zero) {
#if defined(DEBUG)
                        cout << "addBase end work.zero" << endl;
#endif
                        return;
                    }
                }
                int pivot = calc_1pos(work.next);
#if defined(DEBUG)
                cout << "addBase pivot = " << dec << pivot << endl;
#endif
                if (pivot >= size) {
                    cout << "pivot > size error pivot = " << dec << pivot
                         << " size = " << dec << size << endl;
                    throw "pivot > size error";
                }
                if (isZero(bases[pivot].next)) {
                    add_state(bases[pivot], work);
#if defined(DEBUG)
                    cout << "addBase end next == 0" << endl;
#endif
                    return;
                }
                add_state(work, bases[pivot]);
            }
#if defined(DEBUG)
            cout << "addBase end" << endl;
#endif
        }

        /**
         *\japanese
         * qによって殲滅される部分空間の基底からパリティチェック用の定数を求める
         * @param[in] base 基底
         * @param[in] size 基底の数
         * @return パリティチェックベクトル
         *\endjapanese
         *
         *\english
         * Search parameters for state transition function and
         * parameters for tempering.
         * @param[in] base basis
         * @param[in] size size of basis
         * @return a parity check vector
         *\endenglish
         */
        U search_parity_check_vector(internal_state base[], int size) {
#if defined(DEBUG)
            cout << "search_parity_check_vector start" << endl;
#endif
            NTL::mat_GF2 mx;
            NTL::mat_GF2 my;

            mx.SetDims(word_width, size);
            U mask;
            int pos = bit_size<U>() - 1;
            setZero(mask);
            setBitOfPos(&mask, pos, 1);
//            mask = ~mask;
//            mask = mask ^ (mask >> 1);
            for (int i = 0; i < word_width; i++) {
                int cnt = 0;
                for (int j = 0; j < word_width; j++) {
                    if (isZero(base[j].next)) {
                        continue;
                    }
                    if (!isZero(mask & base[j].next)) {
                        mx.put(i, cnt, 1);
                    } else {
                        mx.put(i, cnt, 0);
                    }
                    cnt++;
                }
                pos--;
                setZero(mask);
                setBitOfPos(&mask, pos, 1);
//                mask = mask >> 1;
            }
            kernel(my, mx);
            if (my.NumRows() == 0) {
                cout << "parity vector can't find" << endl;
                throw "parity vector can't find";
            }
#if defined(DEBUG)
            cout << "dim mx = "<< mx.NumRows() << endl;
            cout << "-----" << endl;
            for (int i = 0; i < mx.NumRows(); i++) {
                for (int j = 0; j < mx.NumCols(); j++) {
                    cout << mx.get(i, j);
                }
                cout << endl;
            }
            cout << "dim kernel = "<< my.NumRows() << endl;
            cout << "-----" << endl;
            for (int i = 0; i < my.NumRows(); i++) {
                for (int j = 0; j < my.NumCols(); j++) {
                    cout << my.get(i, j);
                }
                cout << endl;
            }
            cout << "search_parity_check_vector end" << endl;
#endif
            return fromGF2Vec<U>(my[0]);
        }

    }; // End of class AlgorithmSearchParity
} // End of name space MTToolBox
#endif // MTTOOLBOX_ALGORITHM_CALCULATE_PARITY_HPP
