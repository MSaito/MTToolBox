#pragma once
#ifndef CALC_PARITY_H
#define CALC_PARITY_H
#include <NTL/GF2X.h>
#include "sfmtsearch.hpp"

MTToolBox::w128_t calc_parity(const MTToolBox::sfmt& dsfmt,
                              const NTL::GF2X& irreducible);
#endif
