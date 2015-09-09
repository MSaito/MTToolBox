#include "printBinary.h"

using namespace NTL;
using namespace std;

void printBinary(FILE *fp, const GF2X& poly)
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
    fflush(fp);
}
