/**
 * @file tinymt32dc.cpp
 *
 * @brief The main function of parameter generator of 32-bit Tiny
 * Mersenne Twister.
 *
 * The functions in this file are simple. They parse the command line
 * options and call all_in_one function which does almost all things.
 * Users can change this file so that it fits to their applications
 * and OS.
 *
 * @author Mutsuo Saito (Hiroshima University)
 * @author Makoto Matsumoto (The University of Tokyo)
 *
 * Copyright (C) 2011 Mutsuo Saito, Makoto Matsumoto,
 * Hiroshima University and The University of Tokyo.
 * All rights reserved.
 *
 * The 3-clause BSD License is applied to this software, see
 * LICENSE.txt
 */
#include <stdint.h>
#include <inttypes.h>
#include <time.h>
#include <unistd.h>
#include <string>
#include <sstream>
#include <fstream>
#include <MTToolBox/AlgorithmRecursionAndTempering.hpp>
#include <MTToolBox/Sequential.hpp>
#include "tinymt32search.hpp"
#include "parse_opt.hpp"

using namespace std;
using namespace MTToolBox;
using namespace tinymt;
static const uint32_t sequence_max = 0x7fffffff;

int search(tinymt_options& opt, int count);
const std::string toString(GF2X& poly);

/**
 * parse command line option, and search parameters
 * @param argc number of arguments
 * @param argv value of arguments
 * @return 0 if this ends normally
 */
int main(int argc, char** argv) {
    tinymt_options opt;
    bool parse = parse_opt(opt, argc, argv, sequence_max);
    if (!parse) {
        return -1;
    }
    try {
        return search(opt, opt.count);
    } catch (underflow_error &e) {
        return 0;
    }
}

/**
 * search parameters using all_in_one function in the file search_all.hpp
 * @param opt command line options
 * @param count number of parameters user requested
 * @return 0 if this ends normally
 */
int search(tinymt_options& opt, int count) {
    Sequential<uint32_t> sq(0, opt.start);
    tinymt32 g(opt.uid);

    if (opt.verbose) {
        time_t t = time(NULL);
        cout << "search start at " << ctime(&t);
        cout << "id:" << dec << opt.uid << endl;
    }
    AlgorithmRecursionAndTempering<uint32_t> all(sq);
    st32 st;
    stlsb32 stlsb;
    int i = 0;
    while (i < count || opt.all) {
        if (all.search(g, st, stlsb, opt.verbose)) {
            int delta = all.getDelta();
            if (delta > opt.max_delta) {
                continue;
            }
            //g = all.get_rand();
            //tinymt32_param param = g.get_param();
            int weight = all.getWeight();
            GF2X poly = all.getCharacteristicPolynomial();
            if (i == 0) {
                cout << "# characteristic, "
                     << g.getHeaderString()
                     << ", weight, delta"
                     << endl;
            }
            cout << toString(poly) << ","
                 << g.getParamString()
                 << dec << weight << "," << delta
                 << endl;
#if 0
            output_params<uint32_t, tinymt32_param>(poly, weight,
                                                    delta, param,
                                                    opt, i == 0);
#endif
            i++;
        } else {
            cout << "search failed" << endl;
            break;
        }
    }
    if (opt.verbose) {
        time_t t = time(NULL);
        cout << "search end at " << ctime(&t) << endl;
    }
    return 0;
}

const std::string toString(GF2X& poly) {
    uint64_t p1;
    uint64_t p2;
    uint64_t p = 0;
    uint64_t mask = 1;
    for(int i = 0; i < 64; i++) {
        if (IsOne(coeff(poly, i))) {
            p |= mask;
        }
        mask <<= 1;
    }
    p2 = p;
    mask = 1;
    p = 0;
    for(int i = 0; i < 64; i++) {
        if (IsOne(coeff(poly, i + 64))) {
            p |= mask;
        }
        mask <<= 1;
    }
    p1 = p;
    stringstream ss;
    ss << hex << setfill('0') << setw(16) << p1
       << hex << setfill('0') << setw(16) << p2;
    return ss.str();
}
