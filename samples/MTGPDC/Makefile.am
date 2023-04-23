noinst_PROGRAMS = MTGP32dc
check_PROGRAMS = test
MTGP32dc_SOURCES = MTGP32dc.cpp MTGP32search.hpp parse_opt.cpp parse_opt.hpp \
mtgp_param.h
test_SOURCES = test_linearity.cpp

DEPENDENCIES = ../../lib/libMTToolBox.la
LDADD = ../../lib/libMTToolBox.la

#MTGP32dc_DEPENDENCIES = ../../lib/libMTToolBox.a
#test_DEPENDENCIES = ../../lib/libMTToolBox.a
AM_LDFLAGS = -L../../lib
#MTGP32dc_LDFLAGS = -L../../lib
#test_LDFLAGS = -L../../lib
#LDADD = -lMTToolBox
#MTGP32dc_LDADD = -lMTToolBox -lntl
#test_LDADD = -lMTToolBox -lntl

warning = -Wall -Wextra -Wsign-compare -Wconversion
stdcmacros = -D__STDC_CONSTANT_MACROS -D__STDC_FORMAT_MACROS
AM_CXXFLAGS = $(warning) $(stdcmacros) -O2 -I../../include
