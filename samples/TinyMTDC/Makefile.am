noinst_PROGRAMS = tinymt32dc
tinymt32dc_SOURCES = parse_opt.cpp parse_opt.hpp tinymt32dc.cpp \
	tinymt32search.hpp

DEPENDENCIES = ../../lib/libMTToolBox.a
AM_LDFLAGS = -L../../lib
LDADD = -lMTToolBox

EXTRA_DIST = README.jp.txt

AM_CXXFLAGS = -Wall -O2 -I../../include -Wextra -D__STDC_CONSTANT_MACROS \
              -D__STDC_FORMAT_MACROS
