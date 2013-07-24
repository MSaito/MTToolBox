
#ifndef ABSTRACT_SEARCHER_HPP
#define ABSTRACT_SEARCHER_HPP
#include <stdint.h>
#include <inttypes.h>
#include <MTToolBox/abstract_generator.hpp>

/* C++的に正しくないプログラム
   わかりやすくて低速なものを後で高速化することはできるが、
   わかりにくくて高速なものを後でわかりやすくすることは難しい。*/

namespace MTTools {
    class U32Searcher : public U32Generator {
    public:
        virtual void add(U32Searcher that) = 0;
        virtual void next() = 0;
        virtual void setZero() = 0;
        virtual boole isZero() = 0;
    };
    class I32Searcher : public I32Generator {
    public:
        virtual int32_t generate() = 0;
        virtual void seed(int32_t value) = 0;
    };
    class U64Searcher : public U64Generator {
    public:
        virtual uint64_t generate() = 0;
        virtual void seed(uint64_t value) = 0;
    };
    class I64Searcher : public I64Generator {
    public:
        virtual int64_t generate() = 0;
        virtual void seed(int64_t value) = 0;
    };
    class U128Searcher : public U128Generator {
        virtual uint128_t generate() = 0;
        virtual void seed(uint128_t value) = 0;
    };
}
#endif //ABSTRACT_SEARCHER_HPP
