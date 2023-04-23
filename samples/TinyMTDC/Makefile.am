noinst_PROGRAMS = tinymt32dc
tinymt32dc_SOURCES = parse_opt.cpp parse_opt.hpp tinymt32dc.cpp \
	tinymt32search.hpp

DEPENDENCIES = ../../lib/libMTToolBox.la
#LDADD = -lMTToolBox
LDADD = ../../lib/libMTToolBox.la

#AM_LDFLAGS = -L../../lib
EXTRA_DIST = README.jp.txt

warning = -Wall -Wextra -Wsign-compare -Wconversion
stdcmacros = -D__STDC_CONSTANT_MACROS -D__STDC_FORMAT_MACROS
AM_CXXFLAGS = $(warning) $(stdcmacros) -O2 -I../../include

#AM_CXXFLAGS = -Wall -O2 -I../../include -Wextra -D__STDC_CONSTANT_MACROS \
#              -D__STDC_FORMAT_MACROS
