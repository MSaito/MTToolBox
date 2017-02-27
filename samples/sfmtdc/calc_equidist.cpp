#include "sfmtsearch.hpp"
#include "AlgorithmSIMDEquidistribution.hpp"
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
    sfmt_param params;
};

static bool parse_opt(options& opt, int argc, char **argv);
static void output_help(string& pgm);
//static void printBinary(FILE *fp, GF2X& poly);

int main(int argc, char * argv[])
{
    options opt;
    if (!parse_opt(opt, argc, argv)) {
        return -1;
    }
    sfmt sf(opt.params);
    //cout << sf.getParamString() << endl;
    //w128_t wseed = {{0x1234, 0x5678, 0x9abc, 0xdef0}};
    //w128_t wseed = {{0x12345, 0x56789, 0x9abcd, 0xdef01}};
    w128_t wseed;
    wseed.u64[0] = opt.seed;
    sf.seed(wseed);
    if (!anni(sf)) {
        return -1;
    }
    int delta128 = 0;
    int delta32 = 0;
    int delta64 = 0;
    bool lsb = false;
    const char * lsb_str = "";
    if (opt.reverse) {
        sf.set_reverse_bit();
        lsb = true;
        cout << "Equidistribution from LSB" << endl;
        lsb_str = " from LSB";
    }
    AlgorithmEquidistribution<w128_t> re(sf, 128, opt.params.mexp);
    int veq[128];
    delta128 = re.get_all_equidist(veq);
    SIMDInfo info;
    info.bitSize = 128;
    int veq32[32];
    info.bitMode = 32;
    info.elementNo = 4;
    sf.reset_reverse_bit();
    delta32 = calc_SIMD_equidistribution<w128_t, sfmt>(sf, veq32, 32, info,
                                                       opt.params.mexp, lsb);
    int veq64[64];
    info.bitMode = 64;
    info.elementNo = 2;
    delta64 = calc_SIMD_equidistribution<w128_t, sfmt>(sf, veq64, 64, info,
                                                       opt.params.mexp, lsb);
    cout << sf.getParamString();
    cout << dec << delta32 << "," << delta64 << ","
         << delta128 << endl;
    if (opt.verbose) {
        cout << "32bit dimension of equidistribution at v-bit accuracy k(v)"
             << lsb_str << endl;
        for (int j = 0; j < 32; j++) {
            cout << "k(" << dec << (j + 1) << ") = " << dec << veq32[j];
            cout << "\td(" << dec << (j + 1) << ") = " << dec
                 << (opt.params.mexp / (j + 1) - veq32[j]) << endl;
        }
        cout << "64bit dimension of equidistribution at v-bit accuracy k(v)"
             << lsb_str << endl;
        for (int j = 0; j < 64; j++) {
            cout << "k(" << dec << (j + 1) << ") = " << dec << veq64[j];
            cout << "\td(" << dec << (j + 1) << ") = " << dec
                 << (opt.params.mexp / (j + 1) - veq64[j]) << endl;
        }
        cout << "128bit dimension of equidistribution at v-bit accuracy k(v)"
             << lsb_str << endl;
        for (int j = 0; j < 128; j++) {
            cout << "k(" << dec << (j + 1) << ") = " << dec << veq[j];
            cout << "\td(" << dec << (j + 1) << ") = " << dec
                 << (opt.params.mexp / (j + 1) - veq[j]) << endl;
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
    opt.reverse = false;
    int c;
    bool error = false;
    string pgm = argv[0];
    static struct option longopts[] = {
        {"verbose", no_argument, NULL, 'v'},
        {"seed", required_argument, NULL, 's'},
        {"reverse", no_argument, NULL, 'r'},
        {NULL, 0, NULL, 0}};
    errno = 0;
    for (;;) {
        c = getopt_long(argc, argv, "vrs:", longopts, NULL);
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
        case 'r':
            opt.reverse = true;
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
        opt.params.sl2 = strtoul(para, &para, 10);
        para++;
        opt.params.sr1 = strtoul(para, &para, 10);
        para++;
        opt.params.sr2 = strtoul(para, &para, 10);
        para++;
        opt.params.msk1 = strtoul(para, &para, 16);
        para++;
        opt.params.msk2 = strtoul(para, &para, 16);
        para++;
        opt.params.msk3 = strtoul(para, &para, 16);
        para++;
        opt.params.msk4 = strtoul(para, &para, 16);
        para++;
        opt.params.parity1 = strtoul(para, &para, 16);
        para++;
        opt.params.parity2 = strtoul(para, &para, 16);
        para++;
        opt.params.parity3 = strtoul(para, &para, 16);
        para++;
        opt.params.parity4 = strtoul(para, &para, 16);
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
         << " mexp,pos1,sl1,sl2,sr1,sr2,msk1,msk2,msk3,msk4,"
         << "parity1,parity2,parity3,parity4"
         << endl;
    static string help_string1 = "\n"
"--verbose, -v        Verbose mode. Output detailed information.\n"
"--reverse, -r        Reverse mode. Calculate equidistribution from LSB.\n"
        ;
    cerr << help_string1 << endl;
}
