# -*- Makefile -*-
WARN = -Wall -Wextra
#UNROLL = -funroll-loops
OPTI = -O2 #-msse4.2 -mssse3 -msse2 $(UNROLL)
#DEBUG = -DDEBUG -g -ggdb #-O0
#DEBUG = -DDEBUG
#SANITY = -fsanitize=address -fno-omit-frame-pointer -g
CXX = g++
#CXX = clang++
#OBJS = period.o AlgorithmPrimitivity.o
STD = -std=c++11
INCLUDE = -I$(HOME)/git/MTToolbox/include
LIB = -lntl -lgmp -lgf2x -lMTToolBox
CXXFLAGS = -D__STDC_CONSTANT_MACROS -D__STDC_FORMAT_MACROS \
$(OPTI) $(WARN) $(STD) $(DEBUG) $(INCLUDE) $(SANITY) $(LIB)
OBJS = cpu_feature.o aligned_alloc.o

all: test

calc_equidist: calc_equidist.cpp sfmtsearch.hpp \
	AlgorithmSIMDEquidistribution.hpp Annihilate.h Annihilate.o
	$(CXX) $(CXXFLAGS) -o $@ calc_equidist.cpp Annihilate.o

test_linearity: test_linearity.cpp dSFMTsearch.hpp
	$(CXX) $(CXXFLAGS) -o $@ test_linearity.cpp

calc_period: calc_period.cpp dSFMTsearch.hpp Annihilate.o
	$(CXX) $(CXXFLAGS) -o $@ calc_period.cpp Annihilate.o

calc_fixpoint: calc_fixpoint.cpp calc_fixpoint.h dSFMTsearch.hpp \
	Annihilate.o printBinary.o
	$(CXX) $(CXXFLAGS) -DMAIN -o $@ calc_fixpoint.cpp Annihilate.o \
	printBinary.o

make_old_format: make_old_format.cpp dSFMTsearch.hpp printBinary.o Annihilate.o
	$(CXX) $(CXXFLAGS) -o $@ make_old_format.cpp Annihilate.o printBinary.o


sfmtdc: sfmtdc.cpp sfmtsearch.hpp AlgorithmSIMDEquidistribution.hpp \
	Annihilate.o Annihilate.h
	$(CXX) $(CXXFLAGS) -o $@ sfmtdc.cpp Annihilate.o

sfmtdcfixed: sfmtdcfixed.cpp sfmtsearch.hpp AlgorithmSIMDEquidistribution.hpp \
	Annihilate.o Annihilate.h
	$(CXX) $(CXXFLAGS) -o $@ sfmtdcfixed.cpp Annihilate.o

Annihilate.o: dSFMTsearch.hpp Annihilate.cpp Annihilate.h
	$(CXX) $(CXXFLAGS) -c Annihilate.cpp

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -rf *.o *~ *.dSYM