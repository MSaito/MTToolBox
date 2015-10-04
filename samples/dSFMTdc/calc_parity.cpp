#include "dSFMTsearch.hpp"
#include "Annihilate.h"
#include <errno.h>
#include <stdlib.h>
#include <getopt.h>
#include <errno.h>
#include <NTL/GF2X.h>
#include <NTL/GF2XFactoring.h>
#include <NTL/vector.h>
#include <MTToolBox/period.hpp>
#include <MTToolBox/AlgorithmReducibleRecursionSearch.hpp>
#include <MTToolBox/AlgorithmCalculateParity.hpp>
#include "calc_parity.h"

using namespace MTToolBox;
using namespace NTL;
using namespace std;

w128_t calc_parity(dSFMT& dsfmt, const GF2X& irreducible)
{
    AlgorithmCalculateParity<w128_t, dSFMT> cp;
    cp.searchParity(dsfmt, irreducible);
    return dsfmt.getParityValue();
}

