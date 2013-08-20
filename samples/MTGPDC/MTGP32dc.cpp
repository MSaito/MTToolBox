/**
 * @file mtgp32dc.cpp
 *
 * @brief search irreducible polynomial and temper the output for 32bit
 * mtgp.
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
#ifndef __STDC_CONSTANT_MACROS
#define __STDC_CONSTANT_MACROS 1
#endif

#include <stdint.h>
#include <inttypes.h>
#include <time.h>
//#include <unistd.h>
#include <string>
#include <sstream>
#include <fstream>
#include "all_in_one.hpp"
#include "mtgp32.hpp"
#include "parse_opt.h"
#include "output.hpp"

using namespace std;
using namespace MTToolBox;
using namespace mtgp;

typedef linear_generator<uint32_t, mtgp32> lg32;
typedef search_temper<lg32, uint32_t, 32, 4, 23, 5> st32;
typedef search_temper<lg32, uint32_t, 32, 4, 9, 5, true> stlsb32;

int indexed_search(mtgp_options& opt, bool first);

int main(int argc, char** argv) {
    mtgp_options opt;
    bool first = true;
    bool parse = parse_opt(opt, argc, argv);
    if (!parse) {
	return -1;
    }
    while (opt.count > 0) {
	if (indexed_search(opt, first)) {
	    opt.id += 1;
	    opt.count -= 1;
	} else {
	    return -1;
	}
	first = false;
    }
    return 0;
}

int indexed_search(mtgp_options& opt, bool first) {
    mtgp32 mtgp(opt.mexp, opt.id);
    shared_ptr<lg32> g(new lg32(mtgp));

    if (opt.verbose) {
	time_t t = time(NULL);
	cout << "search start at " << ctime(&t);
	cout << "mexp:" << dec << opt.mexp << endl;
	cout << "id:" << dec << opt.id << endl;
    }
    all_in_one<uint32_t, lg32, st32, stlsb32> all;
    st32 st;
    stlsb32 stlsb;
    if (all.search(*g, &st, &stlsb, opt.seed, opt.verbose)) {
	g = all.get_rand();
	shared_ptr<mtgp32> rand = g->get_random();
	mtgp_param<uint32_t> param = rand->get_param();
	int weight = all.get_weight();
	int delta = all.get_delta();
	string sha1 = all.get_sha1();
	param.set_sha1(sha1);
	output_mtgp_params<uint32_t>(weight, delta, param, opt, first);
	return 1;
    } else {
	cout << "search failed" << endl;
	return 0;
    }
}
