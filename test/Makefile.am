check_PROGRAMS = test_mt test_period test_search test_equidist \
	test_temper test_search_temper test_linearity test_primefactors
TESTS = test_mt test_period test_search test_equidist \
	test_temper test_search_temper test_linearity test_primefactors
test_mt_SOURCES = test_mt.cpp tinymt32.c tinymt32.h tinymt64.c tinymt64.h
test_period_SOURCES = test_period.cpp tinymt32.c tinymt32.h tinymt64.c \
	tinymt64.h
test_search_SOURCES = test_search.cpp test_generator.hpp tinymt32.c \
	tinymt32.h tinymt64.c tinymt64.h
test_equidist_SOURCES = test_equidist.cpp tinymt32.c tinymt32.h \
	tinymt64.c tinymt64.h
test_temper_SOURCES = test_temper.cpp tinymt32.c tinymt32.h tinymt64.c \
	tinymt64.h
test_search_temper_SOURCES = test_search_temper.cpp tinymt32.c tinymt32.h \
	tinymt64.c tinymt64.h
test_linearity_SOURCES = test_linearity.cpp tinymt32.c tinymt32.h tinymt64.c \
	tinymt64.h
test_primefactors_SOURCES = test_primefactors.cpp tinymt32.c tinymt32.h \
	tinymt64.c tinymt64.h

#test1_SOURCES = test_period.cpp test_search.cpp test_equidist.cpp
#test2_SOURCES = test_temper.cpp test_search_temper.cpp \
#	test_temper_searcher.hpp test_linearity.cpp test_mt.cpp \
#	test_primefactors.cpp test_generator.hpp
EXTRA_DIST = tinymt32.c tinymt32.h tinymt64.c tinymt64.h

test_mt_DEPENDENCIES = ../lib/libMTToolBox.a
test_period_DEPENDENCIES = ../lib/libMTToolBox.a
test_search_DEPENDENCIES = ../lib/libMTToolBox.a
test_equidist_DEPENDENCIES = ../lib/libMTToolBox.a
test_temper_DEPENDENCIES = ../lib/libMTToolBox.a
test_search_temper_DEPENDENCIES = ../lib/libMTToolBox.a
test_linearity_DEPENDENCIES = ../lib/libMTToolBox.a
test_primefactors_DEPENDENCIES = ../lib/libMTToolBox.a

#test1_DEPENDENCIES = ../lib/libMTToolBox.a
#test2_DEPENDENCIES = ../lib/libMTToolBox.a

test_mt_LDFLAGS =  -L../lib
test_period_LDFLAGS =  -L../lib
test_search_LDFLAGS =  -L../lib
test_equidist_LDFLAGS =  -L../lib
test_temper_LDFLAGS =  -L../lib
test_search_temper_LDFLAGS =  -L../lib
test_linearity_LDFLAGS =  -L../lib
test_primefactors_LDFLAGS =  -L../lib

#test1_LDFLAGS = -L../lib
#test2_LDFLAGS = -L../lib

test_mt_LDADD = -lMTToolBox -lntl
test_period_LDADD = -lMTToolBox -lntl
test_search_LDADD = -lMTToolBox -lntl
test_equidist_LDADD = -lMTToolBox -lntl
test_temper_LDADD = -lMTToolBox -lntl
test_search_temper_LDADD = -lMTToolBox -lntl
test_linearity_LDADD = -lMTToolBox -lntl
test_primefactors_LDADD = -lMTToolBox -lntl

#test1_LDADD = -lMTToolBox -lntl -lUnitTest++
#test2_LDADD = -lMTToolBox -lntl -lUnitTest++

AM_CXXFLAGS = -Wall -I../include -g -Wextra -D__STDC_CONSTANT_MACROS \
              -D__STDC_FORMAT_MACROS
AM_CFLAGS = -Wall -I../include -g -Wextra -Wmissing-prototypes -std=c99
