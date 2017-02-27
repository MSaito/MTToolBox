noinst_PROGRAMS = sfmtdc sfmtdcfixed test_linearity test_period calc_equidist

sfmtdc_SOURCES = sfmtsearch.hpp sfmtdc.cpp Annihilate.h Annihilate.cpp \
AlgorithmSIMDEquidistribution.hpp w128.hpp
sfmtdcfixed_SOURCES = sfmtsearch.hpp sfmtdcfixed.cpp Annihilate.h \
Annihilate.cpp AlgorithmSIMDEquidistribution.hpp w128.hpp
calc_equidist_SOURCES = sfmtsearch.hpp calc_equidist.cpp Annihilate.cpp \
AlgorithmSIMDEquidistribution.hpp
test_linearity_SOURCES = sfmtsearch.hpp test_linearity.cpp Annihilate.cpp
test_period_SOURCES = sfmtsearch.hpp test_period.cpp Annihilate.cpp

DEPENDENCIES = ../../lib/libMTToolBox.a
LDADD = -lMTToolBox
AM_LDFLAGS = -L../../lib
AM_CXXFLAGS = -Wall -O2 -I../../include -Wextra -D__STDC_CONSTANT_MACROS \
              -D__STDC_FORMAT_MACROS
EXTRA_DIST = README.txt README.utf8.jp.txt
