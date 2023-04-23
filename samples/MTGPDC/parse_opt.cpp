/**
 * @file parse_opt.cpp
 *
 * @brief parse command line option
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
#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>
#include <errno.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "parse_opt.hpp"

using namespace std;

bool parse_opt(mtgp_options& opt, int argc, char **argv) {
    int c;
    bool error = false;
    string pgm = argv[0];
    static struct option longopts[] = {
        {"seed", required_argument, NULL, 's'},
        {"verbose", no_argument, NULL, 'v'},
        {"output-file", required_argument, NULL, 'f'},
        {"count", required_argument, NULL, 'c'},
        {NULL, 0, NULL, 0}};
    opt.verbose = false;
    opt.count = 1;
    for (;;) {
        c = getopt_long(argc, argv, "vs:f:c:", longopts, NULL);
        if (error) {
            break;
        }
        if (c == -1) {
            break;
        }
        switch (c) {
        case 'v':
            opt.verbose = true;
            break;
        case 's':
            opt.seed = optarg;
            break;
        case 'f':
            opt.filename = optarg;
            break;
        case 'c':
            opt.count = strtoll(optarg, NULL, 10);
            if (errno) {
                error = true;
                cerr << "count must be number" << endl;
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
    if (argc < 2) {
        error = true;
    } else {
        opt.mexp = static_cast<int>(strtol(argv[0], NULL, 10));
        if (errno != 0) {
            error = true;
            cerr << "mexp must be number" << endl;
        } else {
            if (opt.mexp != 3217
                && opt.mexp != 4423
                && opt.mexp != 11213
                && opt.mexp != 23209
                && opt.mexp != 44497
                && opt.mexp != 110503
                && opt.mexp != 216091) {
                error = true;
                cerr << "mexp is 3217, 4423, 11213, 23209, 44497,"
                    " 110503 or 216091" << endl;
            }
        }
        long long id = strtoll(argv[1], NULL, 10);
        opt.id = static_cast<uint32_t>(id);
        if (errno != 0) {
            error = true;
            cerr << "id must be a number between 0 and 2^32-1" << endl;
        }
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
        cerr << pgm
             << " [-v] [-c count] [-s seed_string] [-f outputfile] mexp id"
             << endl;
        cerr << "mexp                   mersenne exponent of search" << endl;
        cerr << "                       the generator will have the period "
             << "of 2^{mexp}-1." << endl;
        cerr << "index                  index of generator. the parameters searched"
             << endl;
        cerr << "                       by diffrent indexes will generate diffrent sequences." << endl;
        cerr << "--seed,-s seed_string  seed for paramater search." << endl;
        cerr << "                       asci charactor string is used as seed"
             << endl;
        cerr << "--verbose,-v           verbose mode" << endl;
        cerr << "                       output parameters, calculation time, etc."
             << endl;
        cerr << "--output-file,-f filename parameters are output to this file."
             << endl;
        cerr << "                       without this option, output to "
             << "standard output." << endl;
        cerr << "--count,-c count       output count. increment id and "
             << "repeat output." << endl;
        return false;
    }
    if (opt.seed.length() == 0) {
        stringstream ss;
        ss << dec << getpid() << ":";
        ss << dec << clock();
        opt.seed = ss.str();
    }
    return true;
}
