#bin_PROGRAMS = mt1 mt2
noinst_PROGRAMS = mt1 mt2
mt1_SOURCES = mt1.cpp
mt2_SOURCES = mt2.cpp

mt1_DEPENDENCIES = ../../lib/libMTToolBox.a
mt2_DEPENDENCIES = ../../lib/libMTToolBox.a
mt1_LDFLAGS = -L../../lib
mt2_LDFLAGS = -L../../lib
mt1_LDADD = -lMTToolBox -lntl
mt2_LDADD = -lMTToolBox -lntl

#AM_LINKFLAGS = -L ../../lib

AM_CXXFLAGS = -Wall -I../../include -g -O3 -Wextra -D__STDC_CONSTANT_MACROS \
              -D__STDC_FORMAT_MACROS

