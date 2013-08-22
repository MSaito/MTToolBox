#ifndef PARSE_OPT_H
#define PARSE_OPT_H
/**
 * @file parse_opt.h
 *
 * @brief header file of command line parser.
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
#include <string>

class mtgp_options {
public:
    bool verbose;
    std::string seed;
    std::string filename;
    int mexp;
    uint32_t id;
    long long count;
};

bool parse_opt(mtgp_options& opt, int argc, char **argv);

#endif
