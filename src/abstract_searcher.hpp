
#ifndef ABSTRACT_SEARCHER_HPP
#define ABSTRACT_SEARCHER_HPP
#include <stdint.h>
#include <inttypes.h>

/* C++的に正しくないプログラム
   わかりやすくて低速なものを後で高速化することはできるが、
   わかりにくくて高速なものを後でわかりやすくすることは難しい。*/

namespace MTTools {
    class U32Searcher {
    public:
        virtual uint32_t generate() = 0;
        virtual void seed(uint32_t value) = 0;
    };
    class I32Searcher {
    public:
        virtual int32_t generate() = 0;
        virtual void seed(int32_t value) = 0;
    };
    class U64Searcher {
    public:
        virtual uint64_t generate() = 0;
        virtual void seed(uint64_t value) = 0;
    };
    class I64Searcher {
    public:
        virtual int64_t generate() = 0;
        virtual void seed(int64_t value) = 0;
    };
    class U128Searcher {
        virtual uint128_t generate() = 0;
        virtual void seed(uint128_t value) = 0;
    };
}
#endif //ABSTRACT_SEARCHER_HPP
