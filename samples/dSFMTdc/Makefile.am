noinst_PROGRAMS = dSFMTdc calc_parity calc_equidist calc_fixpoint test_linearity
dSFMTdc_SOURCES = dSFMTsearch.hpp dSFMTdc.cpp Annihilate.h Annihilate.cpp \
AlgorithmDSFMTEquidistribution.hpp w128.hpp
dSFMTdc_DEPENDENCIES = ../../lib/libMTToolBox.a
dSFMTdc_LDFLAGS = -L../../lib
dSFMTdc_LDADD = -lMTToolBox -lntl
calc_equidist_SOURCES = dSFMTsearch.hpp calc_equidist.cpp Annihilate.cpp \
AlgorithmDSFMTEquidistribution.hpp
calc_equidist_DEPENDENCIES = ../../lib/libMTToolBox.a
calc_equidist_LDFLAGS = -L../../lib
calc_equidist_LDADD = -lMTToolBox -lntl
test_linearity_SOURCES = dSFMTsearch.hpp test_linearity.cpp Annihilate.cpp
test_linearity_DEPENDENCIES = ../../lib/libMTToolBox.a
test_linearity_LDFLAGS = -L../../lib
test_linearity_LDADD = -lMTToolBox -lntl
calc_fixpoint_SOURCES = dSFMTsearch.hpp calc_fixpoint.cpp Annihilate.cpp \
printBinary.cpp printBinary.h calc_fixpoint.h
calc_fixpoint_DEPENDENCIES = ../../lib/libMTToolBox.a
calc_fixpoint_LDFLAGS = -L../../lib
calc_fixpoint_LDADD = -lMTToolBox -lntl
calc_parity_SOURCES = dSFMTsearch.hpp calc_parity.cpp Annihilate.cpp \
printBinary.cpp printBinary.h calc_fixpoint.h
calc_parity_DEPENDENCIES = ../../lib/libMTToolBox.a
calc_parity_LDFLAGS = -L../../lib
calc_parity_LDADD = -lMTToolBox -lntl

AM_CXXFLAGS = -Wall -O2 -I../../include -Wextra -D__STDC_CONSTANT_MACROS \
              -D__STDC_FORMAT_MACROS
EXTRA_DIST = README.txt README.utf8.jp.txt
