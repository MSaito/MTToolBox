/* -*- coding:utf-8 -*- */
#ifndef PERIOD_HPP
#define PERIOD_HPP

#include <stdint.h>
#include <NTL/GF2X.h>
#include <NTL/vector.h>
#include <MTToolBox/abstract_generator.hpp>

namespace MTToolBox {
    /**
     * T は U32Generator, I32Generator, U64Generator, I64Generator,
     * または U128Generator のどれか
     *
     */
    template<class T> void
    minpoly(NTL::GF2X& poly, T& generator, int pos = 0)
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

    void
    minpoly(NTL::GF2X& poly, U128Generator& generator, int pos = 0);

    bool isIrreducible(const NTL::GF2X& poly);

    bool isPrime(const NTL::GF2X& poly);

    bool isPrime(const NTL::GF2X& poly, int degree,
                 const NTL::Vec<NTL::ZZ>& prime_factors);
    bool isPrime(const NTL::GF2X& poly,
                 int degree, const char * prime_factors[]);
    bool hasFactorOfDegree(NTL::GF2X& poly, long degree);

}
#endif // PERIOD_HPP
