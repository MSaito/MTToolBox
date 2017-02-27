noinst_PROGRAMS = rmt64dc calc_equidist calc_period
rmt64dc_SOURCES = rmt64.hpp rmt64dc.cpp
calc_equidist_SOURCES = calc_equidist.cpp rmt64.hpp
calc_period_SOURCES = calc_period.cpp rmt64.hpp

DEPENDENCIES = ../../lib/libMTToolBox.a
AM_LDFLAGS = -L../../lib
LDADD = -lMTToolBox

AM_CXXFLAGS = -Wall -O2 -I../../include -Wextra -D__STDC_CONSTANT_MACROS \
              -D__STDC_FORMAT_MACROS
