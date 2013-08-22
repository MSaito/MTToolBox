#ifndef OUTPUT_HPP
#define OUTPUT_HPP 1

/**
 * @file output.hpp
 *
 * @brief format the output of mtgps.
 *
 * @author Mutsuo Saito (Hiroshima University)
 * @author Makoto Matsumoto (Hiroshima University)
 *
 * Copyright (c) 2010 Mutsuo Saito, Makoto Matsumoto and Hiroshima
 * University. All rights reserved.
 * Copyright (c) 2011 Mutsuo Saito, Makoto Matsumoto, Hiroshima
 * University and University of Tokyo. All rights reserved.
 *
 * The new BSD License is applied to this software, see LICENSE.txt
 */
#include <stdint.h>
#include <iostream>
#include <fstream>
#include <MTToolBox/AlgorithmRecursionAndTempering.hpp>
#include <MTToolBox/util.hpp>
#include "MTGP32search.hpp"
#include "parse_opt.hpp"

namespace mtgp {
    template<typename T>
    static void output_params_stream(int weight,
                                     int delta,
                                     mtgp_param<T>& param,
                                     bool first,
                                     std::ostream& ost);
    template<typename T>
    void output_mtgp_params(int weight,
                            int delta,
                            mtgp_param<T>& param,
                            mtgp_options& opt,
                            bool first = true) {
        using namespace MTToolBox;
        using namespace std;

        if (!opt.filename.empty()) {
            // should check in parse param
            ofstream ofs(opt.filename.c_str(), std::ios::out | std::ios::app);
            if (ofs) {
                try {
                    output_params_stream<T>(weight, delta, param, first, ofs);
                    ofs.close();
                } catch (...) {
                    ofs.close();
                }
            }
        } else {
            output_params_stream<T>(weight, delta, param, first, cout);
        }
    }

    template<typename T>
    static void output_params_stream(int weight,
                                     int delta,
                                     mtgp_param<T>& param,
                                     bool first,
                                     std::ostream& ost) {
        using namespace MTToolBox;
        using namespace std;

        int bit_size;
        int width;
        if (sizeof(T) == 4) {
            bit_size = 32;
            width = 8;
        } else {
            bit_size = 64;
            width = 16;
        }
        if (first) {
            ost << "# sha1, mexp, type, id, "
                << "pos, sh1, sh2, tbl_0, tbl_1, tbl_2, tbl_3,"
                << "tmp_0, tmp_1, tmp_2, tmp_3, mask, weight, delta" << endl;
        }
        ost << hex;
        ost << '"';
        for (int i = 0; param.sha1[i]; i++) {
            ost << param.sha1[i];
        }
        ost << '"' << ",";
        ost << dec;
        ost << param.mexp << ",";
        if (bit_size == 32) {
            ost << "uint32_t,";
        } else {
            ost << "uint64_t,";
        }
        ost << param.id << ",";
        ost << param.pos << ",";
        ost << param.sh1 << ",";
        ost << param.sh2 << ",";
        ost << hex;
        for (int i = 0; i < 4; i++) {
            ost << "0x" << setw(width) << setfill('0') << param.tbl[i] << ",";
        }
        for (int i = 0; i < 4; i++) {
            ost << "0x" << setw(width) << setfill('0')
                << param.tmp_tbl[i] << ",";
        }
        ost << "0x" << setw(width) << setfill('0') << param.mask << ",";
        ost << dec;
        ost << weight << ",";
        ost << delta << ",";
        ost << hex;
        for (int i = 0; i < 16; i++) {
            ost << "0x" << setw(width) << setfill('0') << param.p[i] << ",";
        }
        for (int i = 0; i < 16; i++) {
            ost << "0x" << setw(width) << setfill('0') << param.tp[i] << ",";
        }
        if (bit_size == 32) {
            for (int i = 0; i < 16; i++) {
                ost << "0x" << setw(width) << setfill('0')
                    << ((param.tp[i] >> 9) | UINT32_C(0x3f800000))
                    << ",";
            }
        } else {
            for (int i = 0; i < 16; i++) {
                ost << "0x" << setw(width) << setfill('0')
                    << ((param.tp[i] >> 12) | UINT64_C(0x3ff0000000000000))
                    << ",";
            }
        }
        ost << endl;
    }
}
#endif
