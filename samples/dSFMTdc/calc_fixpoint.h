#pragma once
#ifndef CALC_FIXPOINT_H
#define CALC_FIXPOINT_H
#include <NTL/GF2X.h>
#include "dSFMTsearch.hpp"

MTToolBox::w128_t calc_fixpoint(const MTToolBox::dSFMT& dsfmt,
                                const NTL::GF2X& irreducible,
                                const NTL::GF2X& quotient);
#endif
