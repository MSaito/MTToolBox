#pragma once
#ifndef ANNIHILATE_H
#define ANNIHILATE_H

#include "dSFMTsearch.hpp"
void getLCMPoly(NTL::GF2X& lcm, const MTToolBox::dSFMT& sf);
bool anni(MTToolBox::dSFMT& sf);

#endif // ANNIHILATE_H
