noinst_PROGRAMS = mt1 mt2
mt1_SOURCES = mt1.cpp
mt2_SOURCES = mt2.cpp

DEPENDENCIES = ../../lib/libMTToolBox.a
LDADD = -lMTToolBox
AM_LDFLAGS = -L../../lib

AM_CXXFLAGS = -Wall -I../../include -O2 -Wextra -D__STDC_CONSTANT_MACROS \
              -D__STDC_FORMAT_MACROS
