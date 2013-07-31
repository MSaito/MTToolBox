/* -*- coding:utf-8 -*- */
#ifndef MTTOOLBOX_ABSTRACT_SEARCHER_HPP
#define MTTOOLBOX_ABSTRACT_SEARCHER_HPP
#include <stdint.h>
#include <inttypes.h>
#include <MTToolBox/abstract_generator.hpp>

namespace MTToolBox {
    template<class U>
    class Searcher : public Generator<U> {
    public:
        virtual U generate(int outBitLen) = 0;
        virtual void add(Searcher& that) = 0;
        virtual void next() = 0;
        virtual void setZero() = 0;
        virtual bool isZero() = 0;
        virtual void setUpParam() = 0;
        virtual void printParam(std::ostream& out) = 0;
    };
}

#endif // MTTOOLBOX_ABSTRACT_SEARCHER_HPP
