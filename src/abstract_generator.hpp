
#ifndef ABSTRACT_GENERATOR_HPP
#define ABSTRACT_GENERATOR_HPP
#include <stdint.h>
#include <inttypes.h>

/* C++的に正しくないプログラム
   わかりやすくて低速なものを後で高速化することはできるが、
   わかりにくくて高速なものを後でわかりやすくすることは難しい。*/

namespace MTTools {
    class U32Generator {
    public:
        virtual uint32_t generate() = 0;
        virtual void seed(uint32_t value) = 0;
    };
    class I32Generator {
    public:
        virtual int32_t generate() = 0;
        virtual void seed(int32_t value) = 0;
    };
    class U64Generator {
    public:
        virtual uint64_t generate() = 0;
        virtual void seed(uint64_t value) = 0;
    };
    class I64Generator {
    public:
        virtual int64_t generate() = 0;
        virtual void seed(int64_t value) = 0;
    };
    class U128Generator {
        virtual uint128_t generate() = 0;
        virtual void seed(uint128_t value) = 0;
    };
}
#endif //ABSTRACT_GENERATOR_HPP
