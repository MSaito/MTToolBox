noinst_PROGRAMS = xorshift1 xorshift2 xorshift3 xorshift4 xorshift5 \
	xorshift128
xorshift1_SOURCES = xorshift-1.cpp
xorshift2_SOURCES = xorshift-2.cpp
xorshift3_SOURCES = xorshift-3.cpp
xorshift4_SOURCES = xorshift-4.cpp
xorshift5_SOURCES = xorshift-5.cpp
xorshift128_SOURCES = xorshift128_main.c xorshift128.c xorshift128.h

xorshift1_DEPENDENCIES = ../../lib/libMTToolBox.a
xorshift2_DEPENDENCIES = ../../lib/libMTToolBox.a
xorshift3_DEPENDENCIES = ../../lib/libMTToolBox.a
xorshift4_DEPENDENCIES = ../../lib/libMTToolBox.a
xorshift5_DEPENDENCIES = ../../lib/libMTToolBox.a
xorshift128_DEPENDENCIES = ../../lib/libMTToolBox.a

xorshift1_LDFLAGS = -L../../lib
xorshift2_LDFLAGS = -L../../lib
xorshift3_LDFLAGS = -L../../lib
xorshift4_LDFLAGS = -L../../lib
xorshift5_LDFLAGS = -L../../lib
xorshift128_LDFLAGS = -L../../lib

xorshift1_LDADD = -lMTToolBox -lntl
xorshift2_LDADD = -lMTToolBox -lntl
xorshift3_LDADD = -lMTToolBox -lntl
xorshift4_LDADD = -lMTToolBox -lntl
xorshift5_LDADD = -lMTToolBox -lntl
xorshift128_LDADD = -lMTToolBox -lntl

EXTRA_DIST = README README.utf8.ja.txt

AM_CXXFLAGS = -Wall -g -O2 -I../../include -Wextra -D__STDC_CONSTANT_MACROS \
              -D__STDC_FORMAT_MACROS
AM_CFLAGS = -Wall -g -O2 -Wextra -std=c99

