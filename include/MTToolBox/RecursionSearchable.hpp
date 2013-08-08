#ifndef MTTOOLBOX_RECURSION_SEARCHABLE_HPP
#define MTTOOLBOX_RECURSION_SEARCHABLE_HPP

/**
 * @file RecursionSearchable.hpp
 * @brief このクラスは状態遷移パラメータ探索を行えるGF(2)線形疑似乱数生成器のクラスである。
 *
 *
 */
#include <stdint.h>
#include <inttypes.h>
#include <MTToolBox/AbstractGenerator.hpp>

namespace MTToolBox {
    /**
     * @class ResursionSearchable
     * @brief このクラスは状態遷移パラメータ探索を行える
     * GF(2)線形疑似乱数生成器のクラスである。
     *
     * @tparam U 疑似乱数生成器の出力のタイプ
     */
    template<class U>
    class RecursionSearchable : public AbstractGenerator<U> {
    public:
        virtual void setUpParam() = 0;
        virtual void printHeader(std::ostream& out) = 0;
        virtual void printParam(std::ostream& out) = 0;
    };
}

#endif // MTTOOLBOX_RECURSION_SEARCHABLE_HPP
