noinst_PROGRAMS = rmt64dc
rmt64dc_SOURCES = rmt64.hpp rmt64dc.cpp
rmt64dc_DEPENDENCIES = ../../lib/libMTToolBox.a
rmt64dc_LDFLAGS = -L../../lib
rmt64dc_LDADD = -lMTToolBox -lntl

AM_CXXFLAGS = -Wall -g -O3 -I../../include -Wextra -D__STDC_CONSTANT_MACROS \
              -D__STDC_FORMAT_MACROS

