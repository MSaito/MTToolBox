#include "dSFMTsearch.hpp"
#include "AlgorithmDSFMTEquidistribution.hpp"
#include <MTToolBox/AlgorithmEquidistribution.hpp>
#include <MTToolBox/AlgorithmReducibleRecursionSearch.hpp>
#include <MTToolBox/period.hpp>
#include <NTL/GF2X.h>
#include <errno.h>
#include <stdlib.h>
#include <getopt.h>
#include "Annihilate.h"

using namespace MTToolBox;
using namespace std;

class options {
public:
    bool verbose;
    bool reverse;
    uint64_t seed;
    dSFMT_param params;
};

static bool parse_opt(options& opt, int argc, char **argv);
static void output_help(string& pgm);

int main(int argc, char * argv[])
{
    options opt;
    if (!parse_opt(opt, argc, argv)) {
        return -1;
    }
    dSFMT sf(opt.params);
    //cout << sf.getParamString() << endl;
    //w128_t wseed = {{0x1234, 0x5678, 0x9abc, 0xdef0}};
    //w128_t wseed = {{0x12345, 0x56789, 0x9abcd, 0xdef01}};
    w128_t wseed;
    wseed.u64[0] = opt.seed;
    sf.seed(wseed);
#if defined(DEBUG)
    sf.d_p();
    cout << endl;
    for (int w = 2; w >= 1; w--) {
        cout << "weight = " << dec << w << endl;
        sf.setWeightMode(w);
        for (int i = 0; i < 10; i++) {
            w128_t tmp = sf.generate();
            cout << hex << setw(16) << setfill('0') << tmp.u64[0] << " ";
            cout << hex << setw(16) << setfill('0') << tmp.u64[1] << endl;
        }
        cout << endl;
        sf.seed(wseed);
    }
    sf.d_p();
#endif
    if (!anni(sf)) {
        return -1;
    }
    int delta52 = 0;
    int veq52[52];
    DSFMTInfo info;
    info.bitSize = 128; // IMPORTANT
    info.elementNo = 2;
    delta52 = calc_dSFMT_equidistribution<w128_t, dSFMT>(sf, veq52, 52, info,
                                                         opt.params.mexp);
    cout << sf.getParamString();
    cout << dec << delta52 << endl;
    if (opt.verbose) {
        cout << "52bit dimension of equidistribution at v-bit accuracy k(v)"
             << endl;
        for (int j = 0; j < 52; j++) {
            cout << "k(" << dec << (j + 1) << ") = " << dec << veq52[j];
            cout << "\td(" << dec << (j + 1) << ") = " << dec
                 << (opt.params.mexp / (j + 1) - veq52[j]) << endl;
        }
    }
    return 0;
}

#if 0
static void printBinary(FILE *fp, GF2X& poly)
{
    int i;
    if (deg(poly) < 0) {
        fprintf(fp, "0deg=-1\n");
        return;
    }
    for(i = 0; i <= deg(poly); i++) {
        if(rep(coeff(poly, i)) == 1) {
            fprintf(fp, "1");
        } else {
            fprintf(fp, "0");
        }
        if ((i % 32) == 31) {
            fprintf(fp, "\n");
        }
    }
    fprintf(fp, "deg=%ld\n", deg(poly));
}
#endif
/**
 * command line option parser
 * @param opt a structure to keep the result of parsing
 * @param argc number of command line arguments
 * @param argv command line arguments
 * @param start default start value
 * @return command line options have error, or not
 */
static bool parse_opt(options& opt, int argc, char **argv) {
    opt.verbose = false;
    opt.seed = (uint64_t)clock();
    int c;
    bool error = false;
    string pgm = argv[0];
    static struct option longopts[] = {
        {"verbose", no_argument, NULL, 'v'},
        {"seed", required_argument, NULL, 's'},
        {NULL, 0, NULL, 0}};
    errno = 0;
    for (;;) {
        c = getopt_long(argc, argv, "vs:", longopts, NULL);
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
        char * para = argv[0];
        opt.params.mexp = strtoul(para, &para, 10);
        para++;
        opt.params.pos1 = strtoul(para, &para, 10);
        para++;
        opt.params.sl1 = strtoul(para, &para, 10);
        para++;
        opt.params.msk1 = strtoull(para, &para, 16);
        para++;
        opt.params.msk2 = strtoull(para, &para, 16);
        para++;
        opt.params.fix1 = strtoull(para, &para, 16);
        para++;
        opt.params.fix2 = strtoull(para, &para, 16);
        para++;
        opt.params.parity1 = strtoull(para, &para, 16);
        para++;
        opt.params.parity2 = strtoull(para, &para, 16);
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
static void output_help(string& pgm)
{
    cerr << "usage:" << endl;
    cerr << pgm
         << " [-v] [-r]"
         << " mexp,pos1,sl1,msk1,msk2,fix1,fix2,"
         << "parity1,parity2"
         << endl;
    static string help_string1 = "\n"
"--verbose, -v        Verbose mode. Output detailed information.\n"
        ;
    cerr << help_string1 << endl;
}
