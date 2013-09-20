#ifndef MTTOOLBOX_SEQUENTIAL_HPP
#define MTTOOLBOX_SEQUENTIAL_HPP
/**
 * @file Sequential.hpp
 *
 * @brief 順番にカウントダウンする数を生成する
 *
 * @author Mutsuo Saito (Hiroshima University)
 * @author Makoto Matsumoto (Hiroshima University)
 *
 * Copyright (c) 2013 Mutsuo Saito, Makoto Matsumoto and Hiroshima
 * University. All rights reserved.
 *
 * The 3-clause BSD License is applied to this software, see
 * LICENSE.txt
 */

#include <stdexcept>
#include <stdint.h>
#include <inttypes.h>
#include <MTToolBox/AbstractGenerator.hpp>
#include <MTToolBox/util.hpp>

namespace MTToolBox {
    /**
     * @class Sequential
     * @brief カウントダウン生成器
     *
     * @tparam T 出力のタイプ
     */
    template<typename T>
    class Sequential : public AbstractGenerator<T> {
    public:
        Sequential() {
        }
        Sequential(T p_mask) {
            status = reinterpret_cast<T>(-1);
            mask = p_mask;
            error = false;
        }
        Sequential(T p_mask, T seed) {
            status = seed;
            mask = p_mask;
            error = false;
        }
        Sequential(Sequential<T>& src) {
            status = src.status;
            mask = src.mask;
            error = src.error;
        }
        void seed(T value) {
            reseed(value);
        }
        void reseed(T seed) {
            status = seed;
            error = false;
        }
        T generate() {
            return next();
        }
        T next() {
            if (error) {
                throw std::underflow_error("count over zero exception");
            }
            if (status <= 0) {
                error = true;
            }
            T work = status;
            status -= 1;
            return work ^ mask;
        }
        int bitSize() const {
            int r = bit_size<T>();
            return r;
        }
    private:
        T status;
        T mask;
        bool error;
    };
}
//  MTTOOLBOX_SEQUENTIAL_HPP

#endif

