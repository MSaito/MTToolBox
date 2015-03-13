#include <sstream>
#include <NTL/GF2XFactoring.h>
#include <NTL/ZZ.h>
#include <NTL/vector.h>
#include <MTToolBox/period.hpp>
#include <ctype.h>


namespace MTToolBox {
    using namespace std;
    using namespace NTL;

    /**
     * メルセンヌ指数 p のリスト
     */
    static const int32_t mersenne_exponent[] =
    {2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279,
     2203, 2281, 3217, 4253, 4423, 9689, 9941, 11213, 19937, 21701,
     23209, 44497, 86243, 110503, 132049, 216091, 756839, 859433,
     1257787, 1398269, 2976221, 3021377, 6972593, 13466917,
     20996011, 25964951, -1};
#if 0
    /**
     * generator の生成する数列の最小多項式を求める
     * pos はビット位置であり、p = 0 はLSB を指定することになる。
     * generator は出力ではないが、疑似乱数を生成するために、内部状態が
     * 更新される。
     * @param[out] poly 最小多項式
     * @param[in] generator F_2 線形疑似乱数生成器
     * @param[in] pos ビット位置
     */
    void
    minpoly(NTL::GF2X& poly, U128Generator& generator, int pos)
    {
        Vec<GF2> v;
        int size = generator.bitSize();
        v.SetLength(2 * size);
        int index;
        if (pos < 64) {
            index = 0;
        } else {
            index = 1;
            pos = pos - 64;
        }
        for (int i = 0; i < 2 * size; i++) {
            uint128_t u = generator.generate();
            v[i] = (u.u64[index] >> pos) & 1;
        }
        MinPolySeq(poly, v, size);
    }
#endif
    /**
     * 既約判定
     * poly が既約であるか判定する。
     * @param[in] poly
     * @returns true polyが既約が公式の場合
     */
    bool
    isIrreducible(const NTL::GF2X& poly)
    {
        return static_cast<bool>(IterIrredTest(poly));
    }

    /**
     * メルセンヌ指数判定
     * degree がメルセンヌ指数であるか判定する。
     * @param[in] degree
     * @returns true degree がメルセンヌ指数の場合
     */
    bool
    isMexp(uint32_t degree) {
        for (int i = 0; mersenne_exponent[i] > 0; i++) {
            if ((int32_t)degree == mersenne_exponent[i]) {
                return true;
            }
        }
        return false;
    }

    /**
     * 原始性判定
     * poly が原始多項式であるか判定する。
     * この関数では、多項式が既約でかつその次数がメルセンヌ指数
     * の場合に原始多項式であると判定している。
     * この関数によってfalseが返された場合でも、実際には原始多項式
     * であるという可能性がある。次数がメルセンヌ指数と分かっている
     * 場合に使用すること。
     * @param[in] poly
     * @returns true poly が原始多項式の場合
     */
    bool
    isPrime(const NTL::GF2X& poly)
    {
        if (!isIrreducible(poly)) {
            return false;
        }
        uint32_t degree = static_cast<uint32_t>(deg(poly));
        if (isMexp(degree)) {
            return true;
        }
        return false;
    }

    /**
     * 原始性判定
     * poly が原始多項式であるか判定する。この関数は、原始性を正しく
     * 判定するが、2<sup>degree</sup> -1 の素因数分解結果
     * prime_factors を必要とする。prime_factors が正しくない
     * と正しく原始性判定ができない。
     * 2<sup>126</sup>-1 = 3<sup>3</sup> 7<sup>2</sup> 19 43 73 127 337
     * 5419 92737 649657 77158673929 であるが、
     * 素数だけリストにすること。重複度または指数部はリストに入れない。
     * @param[in] poly 判定対象多項式
     * @param[in] degree poly に期待する次数
     * @param[in] prime_factors 素因数分解結果
     */
    bool
    isPrime(const NTL::GF2X& poly,
            int degree, const NTL::Vec<NTL::ZZ>& prime_factors)
    {
        if (deg(poly) != degree) {
            return false;
        }
        if (!isIrreducible(poly)) {
            return false;
        }
        long len = prime_factors.length();
        ZZ period;
        period = 2;
        period = power(period, degree);
        period -= 1;
        for (long i = 0; i < len; i++) {
            ZZ p = prime_factors[i];
            ZZ pow = period / p;
            GF2X x;
            PowerXMod(x, pow, poly);
            if (IsOne(x)) {
                return false;
            }
        }
        return true;
    }

    /**
     * 原始性判定
     * poly が原始多項式であるか判定する。この関数は、原始性を正しく
     * 判定するが、2<sup>degree</sup> -1 の素因数分解結果
     * prime_factors を必要とする。prime_factors が正しくない
     * と正しく原始性判定ができない。
     * 2<sup>126</sup>-1 = 3<sup>3</sup> 7<sup>2</sup> 19 43 73 127 337
     * 5419 92737 649657 77158673929 であるが、
     * 素数だけリストにすること。重複度または指数部はリストに入れない。
     * @param[in] poly 判定対象多項式
     * @param[in] degree poly に期待する次数
     * @param[in] prime_factors 素因数分解結果
     */
    bool
    isPrime(const NTL::GF2X& poly,
            int degree, const char * prime_factors[])
    {
        Vec<ZZ> zz_table;
        long length = 0;
        for (int i = 0; prime_factors[i]; i++) {
            length = i + 1;
        }
        zz_table.SetLength(length);
        for (int i = 0; prime_factors[i] != NULL; i++) {
            ZZ w;
#if defined(ISTREAM)
	    stringstream ss;
            ss << prime_factors[i];
            ss << " ";
            ss >> w;
            zz_table[i] = w;
#else
	    for (const char * p = prime_factors[i]; *p != 0; p++) {
		if (!isdigit(*p)) {
		    break;
		}
		w = w * 10 + (*p) - '0';
	    }
	    zz_table[i] = w;
#endif
        }
        return isPrime(poly, degree, zz_table);
    }

    /**
     * 与えられた多項式が degree 次の既約因子を持つか判定する。
     * poly は常に破壊される。結果がtrue の時、polyには
     * 指定された次数の既約多項式がセットされる。
     * degree は一般にはメルセンヌ指数であり、polyの次数より
     * 少しだけ小さい。degree 次の既約多項式が複数
     * 存在するという可能性は考慮していない。
     */
    bool
    hasFactorOfDegree(NTL::GF2X& poly, long degree)
    {
        static const GF2X t2 = GF2X(2, 1);
        static const GF2X t1 = GF2X(1, 1);
        GF2X t2m;
        GF2X t;
        GF2X alpha;
        int m;

        t2m = t2;
        //DPRINT("degree = %u\n", degree);
        //DPRINTPOLY("poly =", poly);
        if (deg(poly) < degree) {
            return 0;
        }
        t = t1;
        t += t2m;

        for (m = 1; deg(poly) > degree; m++) {
            //DPRINTPOLY("t =", t);
            for(;;) {
                GCD(alpha, poly, t);
                //DPRINTPOLY("alpha =", alpha);
                if (IsOne(alpha)) {
                    break;
                }
                poly /= alpha;
                //DPRINTPOLY("f =", poly);
                if (deg(poly) < degree) {
                    return 0;
                }
            }
            if ((deg(poly) > degree) && (deg(poly) <= degree + m)) {
                //DPRINT("maybe poly is larger m = %d, DEG = %u\n", m,
                //     (unsigned int)deg(poly));
                return 0;
            }
            t2m *= t2m;
            t2m %= poly;
            add(t, t2m, t1);
        }
        if (deg(poly) != degree) {
            return 0;
        }
        return static_cast<bool>(IterIrredTest(poly));
    }

}
