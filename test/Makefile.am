common_files = tinymt32.c tinymt32.h tinymt64.c tinymt64.h
check_PROGRAMS = test_mt test_period test_search test_equidist \
	test_temper test_search_temper test_linearity test_primefactors \
	test_version
TESTS = test_mt test_period test_search test_equidist \
	test_temper test_search_temper test_linearity test_primefactors \
	test_version

test_mt_SOURCES = test_mt.cpp $(common_files)
test_period_SOURCES = test_period.cpp $(common_files)
test_search_SOURCES = test_search.cpp test_generator.hpp $(common_files)
test_equidist_SOURCES = test_equidist.cpp $(common_files)
test_temper_SOURCES = test_temper.cpp $(common_files)
test_search_temper_SOURCES = test_search_temper.cpp $(common_files)
test_linearity_SOURCES = test_linearity.cpp $(common_files)
test_primefactors_SOURCES = test_primefactors.cpp $(common_files)
test_version_SOURCES = test_version.cpp

EXTRA_DIST = $(common_files)

DEPENDENCIES = ../lib/libMTToolBox.la

#LDADD = -lMTToolBox
LDADD = ../lib/libMTToolBox.la

AM_LDFLAGS =  -L../lib -static

warning = -Wall -Wextra -Wsign-compare -Wconversion
stdcmacros = -D__STDC_CONSTANT_MACROS -D__STDC_FORMAT_MACROS
AM_CXXFLAGS = $(warning) $(stdcmacros) -I../include -g
AM_CFLAGS = -Wall -I../include -g -Wextra -Wmissing-prototypes -std=c99
