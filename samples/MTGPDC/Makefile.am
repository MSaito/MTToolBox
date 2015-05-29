noinst_PROGRAMS = MTGP32dc
#noinst_PROGRAMS = MTGP32dc test
MTGP32dc_SOURCES = MTGP32dc.cpp MTGP32search.hpp parse_opt.cpp parse_opt.hpp
#test_SOURCES = test.cpp test_linearity.cpp

MTGP32dc_DEPENDENCIES = ../../lib/libMTToolBox.a
#test_DEPENDENCIES = ../../lib/libMTToolBox.a
MTGP32dc_LDFLAGS = -L../../lib
#test_LDFLAGS = -L../../lib
MTGP32dc_LDADD = -lMTToolBox -lntl
#test_LDADD = -lMTToolBox -lntl

AM_CXXFLAGS = -Wall -I../../include -g -O3 -Wextra -D__STDC_CONSTANT_MACROS \
              -D__STDC_FORMAT_MACROS

EXTRA_DIST = mtgp_param.h

