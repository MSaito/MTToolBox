/* -*- coding:utf-8 -*- */
#ifndef MTTOOLBOX_PERIOD_HPP
#define MTTOOLBOX_PERIOD_HPP

#include <stdint.h>
#include <NTL/GF2X.h>
#include <NTL/vector.h>
#include <MTToolBox/AbstractGenerator.hpp>

namespace MTToolBox {
    /**
     * U は 疑似乱数生成器の出力の型
     *
     */
    template<typename U> void
    minpoly(NTL::GF2X& poly, AbstractGenerator<U>& generator, int pos = 0)
    {
        using namespace std;
        using namespace NTL;

        Vec<GF2> v;
        int size = generator.bitSize();
        v.SetLength(2 * size);
        for (int i = 0; i < 2 * size; i++) {
            v[i] = (generator.generate() >> pos) & 1;
        }
        MinPolySeq(poly, v, size);
    }

    bool isMexp(uint32_t degree);
    bool isIrreducible(const NTL::GF2X& poly);

    bool isPrime(const NTL::GF2X& poly);

    bool isPrime(const NTL::GF2X& poly, int degree,
                 const NTL::Vec<NTL::ZZ>& prime_factors);
    bool isPrime(const NTL::GF2X& poly,
                 int degree, const char * prime_factors[]);
    bool hasFactorOfDegree(NTL::GF2X& poly, long degree);

}
#endif // MTTOOLBOX_PERIOD_HPP
