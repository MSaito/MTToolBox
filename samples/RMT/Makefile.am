noinst_PROGRAMS = rmt64dc calc_equidist calc_period
rmt64dc_SOURCES = rmt64.hpp rmt64dc.cpp
rmt64dc_DEPENDENCIES = ../../lib/libMTToolBox.a
rmt64dc_LDFLAGS = -L../../lib
rmt64dc_LDADD = -lMTToolBox -lntl
calc_equidist_SOURCES = calc_equidist.cpp rmt64.hpp
calc_equidist_DEPENDENCIES = ../../lib/libMTToolBox.a
calc_equidist_LDFLAGS = -L../../lib
calc_equidist_LDADD = -lMTToolBox -lntl
calc_period_SOURCES = calc_period.cpp rmt64.hpp
calc_period_DEPENDENCIES = ../../lib/libMTToolBox.a
calc_period_LDFLAGS = -L../../lib
calc_period_LDADD = -lMTToolBox -lntl

AM_CXXFLAGS = -Wall -O2 -I../../include -Wextra -D__STDC_CONSTANT_MACROS \
              -D__STDC_FORMAT_MACROS

