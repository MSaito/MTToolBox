#ifndef MTTOOLBOX_ABSTRACT_GENERATOR_HPP
#define MTTOOLBOX_ABSTRACT_GENERATOR_HPP
#include <stdint.h>
#include <inttypes.h>
#include <MTToolBox/uint128.hpp>

/* C++的に正しいプログラム
 */

namespace MTToolBox {
    template<class U>
    class Generator {
    public:
        virtual U generate() = 0;
        virtual void seed(U value) = 0;
        virtual int bitSize() const = 0;
    };
}
#endif //MTTOOLBOX_ABSTRACT_GENERATOR_HPP
