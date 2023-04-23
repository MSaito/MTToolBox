noinst_PROGRAMS = mt1 mt2
mt1_SOURCES = mt1.cpp
mt2_SOURCES = mt2.cpp

DEPENDENCIES = ../../lib/libMTToolBox.la

#LDADD = -lMTToolBox
LDADD = ../../lib/libMTToolBox.la

#AM_LDFLAGS = -L../../lib

warning = -Wall -Wextra -Wsign-compare -Wconversion
stdcmacros = -D__STDC_CONSTANT_MACROS -D__STDC_FORMAT_MACROS
AM_CXXFLAGS = $(warning) $(stdcmacros) -O2 -I../../include
#AM_CXXFLAGS = -Wall -Wextra -I../../include -O2 -D__STDC_CONSTANT_MACROS \
#              -D__STDC_FORMAT_MACROS
