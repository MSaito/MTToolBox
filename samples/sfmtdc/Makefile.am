noinst_PROGRAMS = sfmtdc sfmtdcfixed test_linearity test_period calc_equidist
sfmtdc_SOURCES = sfmtsearch.hpp sfmtdc.cpp Annihilate.h Annihilate.cpp \
AlgorithmSIMDEquidistribution.hpp w128.hpp
sfmtdc_DEPENDENCIES = ../../lib/libMTToolBox.a
sfmtdc_LDFLAGS = -L../../lib
sfmtdc_LDADD = -lMTToolBox -lntl
sfmtdcfixed_SOURCES = sfmtsearch.hpp sfmtdcfixed.cpp Annihilate.h \
Annihilate.cpp AlgorithmSIMDEquidistribution.hpp w128.hpp
sfmtdcfixed_DEPENDENCIES = ../../lib/libMTToolBox.a
sfmtdcfixed_LDFLAGS = -L../../lib
sfmtdcfixed_LDADD = -lMTToolBox -lntl
calc_equidist_SOURCES = sfmtsearch.hpp calc_equidist.cpp Annihilate.cpp \
AlgorithmSIMDEquidistribution.hpp
calc_equidist_DEPENDENCIES = ../../lib/libMTToolBox.a
calc_equidist_LDFLAGS = -L../../lib
calc_equidist_LDADD = -lMTToolBox -lntl
test_linearity_SOURCES = sfmtsearch.hpp test_linearity.cpp Annihilate.cpp
test_linearity_DEPENDENCIES = ../../lib/libMTToolBox.a
test_linearity_LDFLAGS = -L../../lib
test_linearity_LDADD = -lMTToolBox -lntl
test_period_SOURCES = sfmtsearch.hpp test_period.cpp Annihilate.cpp
test_period_DEPENDENCIES = ../../lib/libMTToolBox.a
test_period_LDFLAGS = -L../../lib
test_period_LDADD = -lMTToolBox -lntl

AM_CXXFLAGS = -Wall -O2 -I../../include -Wextra -D__STDC_CONSTANT_MACROS \
              -D__STDC_FORMAT_MACROS
EXTRA_DIST = README.txt README.utf8.jp.txt
