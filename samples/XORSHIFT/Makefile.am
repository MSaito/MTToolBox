noinst_PROGRAMS = xorshift1 xorshift2 xorshift3 xorshift4 xorshift5 \
	xorshift128
xorshift1_SOURCES = xorshift-1.cpp
xorshift2_SOURCES = xorshift-2.cpp
xorshift3_SOURCES = xorshift-3.cpp
xorshift4_SOURCES = xorshift-4.cpp
xorshift5_SOURCES = xorshift-5.cpp
xorshift128_SOURCES = xorshift128_main.c xorshift128.c xorshift128.h

DEPENDENCIES = ../../lib/libMTToolBox.la
#LDADD = -lMTToolBox
LDADD = ../../lib/libMTToolBox.la

#AM_LDFLAGS = -L../../lib

EXTRA_DIST = README README.utf8.ja.txt

warning = -Wall -Wextra -Wsign-compare -Wconversion
stdcmacros = -D__STDC_CONSTANT_MACROS -D__STDC_FORMAT_MACROS
AM_CXXFLAGS = $(warning) $(stdcmacros) -O2 -I../../include

AM_CFLAGS = -Wall -g -O2 -Wextra -std=c99
