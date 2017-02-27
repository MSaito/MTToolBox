/**
 * @file sfmtdc.cpp
 *
 * @brief The main function of parameter generator of SFMT.
 *
 * The functions in this file are simple. They parse the command line
 * options and call all_in_one function which does almost all things.
 * Users can change this file so that it fits to their applications
 * and OS.
 *
 * @author Mutsuo Saito (Manieth corp.)
 * @author Makoto Matsumoto (Hiroshima University)
 *
 * Copyright (C) 2015 Mutsuo Saito, Makoto Matsumoto,
 * Hiroshima University.
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
#include <MTToolBox/AlgorithmReducibleRecursionSearch.hpp>
#include <MTToolBox/AlgorithmCalculateParity.hpp>
#include <MTToolBox/AlgorithmEquidistribution.hpp>
#include <MTToolBox/MersenneTwister.hpp>
#include <NTL/GF2X.h>
#include <getopt.h>
#include "sfmtsearch.hpp"
#include "AlgorithmSIMDEquidistribution.hpp"
#include "Annihilate.h"

using namespace std;
using namespace MTToolBox;
using namespace NTL;

class options {
public:
    int mexp;
    bool verbose;
    uint64_t seed;
    std::string filename;
    long count;
};

bool parse_opt(options& opt, int argc, char **argv);


int search(options& opt, int count);

/**
 * parse command line option, and search parameters
 * @param argc number of arguments
 * @param argv value of arguments
 * @return 0 if this ends normally
 */
int main(int argc, char** argv) {
    options opt;
    bool parse = parse_opt(opt, argc, argv);
    if (!parse) {
        return -1;
    }
    return search(opt, opt.count);
}

/**
 * search parameters using all_in_one function in the file search_all.hpp
 * @param opt command line options
 * @param count number of parameters user requested
 * @return 0 if this ends normally
 */
int search(options& opt, int count) {
    MersenneTwister mt(opt.seed);
    sfmt g(opt.mexp);

    cout << "seed = " << dec << opt.seed << endl;
    if (opt.verbose) {
        time_t t = time(NULL);
        cout << "search start at " << ctime(&t);
    }
    AlgorithmReducibleRecursionSearch<w128_t> ars(g, mt);
    int i = 0;
    AlgorithmCalculateParity<w128_t, sfmt> cp;
    cout << "# " << g.getHeaderString() << ", delta32, delta64, delta128"
         << endl;
    while (i < count) {
        if (ars.start(opt.mexp * 100)) {
            GF2X irreducible = ars.getIrreducibleFactor();
            GF2X characteristic = ars.getCharacteristicPolynomial();
            GF2X quotient = characteristic / irreducible;
            //cout << "deg irreducible = " << dec << deg(irreducible) << endl;
            //cout << "deg characteristic = " << dec << deg(characteristic)
            //     << endl;
            //cout << "deg quotient = " << dec << deg(quotient) << endl;
            if (deg(irreducible) != opt.mexp) {
                cout << "error" << endl;
                return -1;
            }
            cp.searchParity(g, irreducible);
            w128_t seed = {{1, 0, 0, 0}};
            g.seed(seed);
            if (!anni(g)) {
                return -1;
            }
            //annihilate<w128_t>(&g, quotient);
            //AlgorithmReducibleEquidistribution<w128_t, sfmt, uint32_t>
            //    re(g, irreducible, 128);
            AlgorithmEquidistribution<w128_t> re(g, 128, opt.mexp);
            int veq[128];
            int delta128 = 0;
            int delta32 = 0;
            int delta64 = 0;
            delta128 = re.get_all_equidist(veq);
            int veq32[32];
            SIMDInfo info;
            info.bitMode = 32;
            info.elementNo = 4;
            info.bitSize = 128;
            delta32
                = calc_SIMD_equidistribution<w128_t, sfmt>(g, veq32, 32, info,
                                                           opt.mexp);
            int veq64[64];
            info.bitMode = 64;
            info.elementNo = 2;
            delta64
                = calc_SIMD_equidistribution<w128_t, sfmt>(g, veq64, 64, info,
                                                           opt.mexp);
            cout << g.getParamString();
            cout << dec << delta32 << "," << delta64 << ","
                 << delta128 << endl;
#if 0
            cout << "32bit veq" << endl;
            for (int j = 0; j < 32; j++) {
                cout << "k(" << dec << (j + 1) << ") = "
                     << dec << veq32[j] << endl;
            }
            cout << "64bit veq" << endl;
            for (int j = 0; j < 64; j++) {
                cout << "k(" << dec << (j + 1) << ") = "
                     << dec << veq64[j] << endl;
            }
            cout << "1282bit veq" << endl;
            for (int j = 0; j < 128; j++) {
                cout << "k(" << dec << (j + 1)<< ") = "
                     << dec << veq[j] << endl;
            }
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

static void output_help(string& pgm);

/**
 * command line option parser
 * @param opt a structure to keep the result of parsing
 * @param argc number of command line arguments
 * @param argv command line arguments
 * @param start default start value
 * @return command line options have error, or not
 */
bool parse_opt(options& opt, int argc, char **argv) {
    opt.verbose = false;
    opt.mexp = 0;
    opt.count = 1;
    opt.seed = (uint64_t)clock();
    opt.filename = "";
    int c;
    bool error = false;
    string pgm = argv[0];
    static struct option longopts[] = {
        {"verbose", no_argument, NULL, 'v'},
        {"file", required_argument, NULL, 'f'},
        {"count", required_argument, NULL, 'c'},
        {"seed", required_argument, NULL, 's'},
        {NULL, 0, NULL, 0}};
    errno = 0;
    for (;;) {
        c = getopt_long(argc, argv, "vs:f:c:", longopts, NULL);
        if (error) {
            break;
        }
        if (c == -1) {
            break;
        }
        switch (c) {
        case 's':
            opt.seed = strtoull(optarg, NULL, 0);
            if (errno) {
                error = true;
                cerr << "seed must be a number" << endl;
            }
            break;
        case 'v':
            opt.verbose = true;
            break;
        case 'f':
            opt.filename = optarg;
            break;
        case 'c':
            opt.count = strtoll(optarg, NULL, 10);
            if (errno) {
                error = true;
                cerr << "count must be a number" << endl;
            }
            break;
        case '?':
        default:
            error = true;
            break;
        }
    }
    argc -= optind;
    argv += optind;
    if (argc < 1) {
        error = true;
    } else {
        long mexp = strtol(argv[0], NULL, 10);
        static const int allowed_mexp[] = {607, 1279,
                                           2281, 3217, 4253, 4423,
                                           9689, 9941, 11213, 19937,
                                           21701, 23209, 44497, -1};
        if (! errno) {
            bool found = false;
            for (int i = 0; allowed_mexp[i] > 0; i++) {
                if (mexp == allowed_mexp[i]) {
                    found = true;
                    break;
                }
            }
            if (! found) {
                error = true;
            }
        }
        if (errno || error){
            error = true;
            cerr << "mexp must be one of ";
            for (int i = 0; allowed_mexp[i] > 0; i++) {
                cerr << dec << allowed_mexp[i] << " ";
            }
            cerr << endl;
        }
        opt.mexp = mexp;
    }
    if (!opt.filename.empty()) {
        ofstream ofs(opt.filename.c_str());
        if (ofs) {
            ofs.close();
        } else {
            error = true;
            cerr << "can't open file:" << opt.filename << endl;
        }
    }
    if (error) {
        output_help(pgm);
        return false;
    }
    return true;
}

/**
 * showing help message
 * @param pgm program name
 */
static void output_help(string& pgm) {
    cerr << "usage:" << endl;
    cerr << pgm
         << " [-s seed] [-v] [-c count]"
         << " [-f outputfile]"
         << " mexp"
         << endl;
    static string help_string1 = "\n"
"--verbose, -v        Verbose mode. Output parameters, calculation time, etc.\n"
"--file, -f filename  Parameters are outputted to this file. without this\n"
"                     option, parameters are outputted to standard output.\n"
"--count, -c count    Output count. The number of parameters to be outputted.\n"
"--seed, -s seed      seed of randomness.\n"
"mexp                 mersenne exponent.\n"
        ;
    cerr << help_string1 << endl;
}
