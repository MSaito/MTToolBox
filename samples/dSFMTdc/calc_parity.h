#pragma once
#ifndef CALC_PARITY_H
#define CALC_PARITY_H
#include <NTL/GF2X.h>
#include "dSFMTsearch.hpp"

MTToolBox::w128_t calc_parity(MTToolBox::dSFMT& dsfmt,
                              const NTL::GF2X& irreducible);
#endif
