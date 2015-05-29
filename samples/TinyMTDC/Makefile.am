noinst_PROGRAMS = tinymt32dc
tinymt32dc_SOURCES = parse_opt.cpp parse_opt.hpp tinymt32dc.cpp \
	tinymt32search.hpp
tinymt32dc_DEPENDENCIES = ../../lib/libMTToolBox.a
tinymt32dc_LDFLAGS = -L../../lib
tinymt32dc_LDADD = -lMTToolBox -lntl

EXTRA_DIST = README.jp.txt

AM_CXXFLAGS = -Wall -g -O3 -I../../include -Wextra -D__STDC_CONSTANT_MACROS \
              -D__STDC_FORMAT_MACROS

