#for GNU make

#DDEBUG = -O0 -g -ggdb -DDEBUG=1

# uncomment next line if you installed NTL with gf2x or gmp
#OPTLIB = -lgf2x
#OPTLIB = -lgmp

INC = -I../include

CC = gcc -Wall -O3 -std=c99 -Wmissing-prototypes -Wextra $(DDEBUG)

CPPOPT = -Wall -Wextra -O3 $(M64) $(INC) \
	$(DDEBUG) -D__STDC_CONSTANT_MACROS -D__STDC_FORMAT_MACROS
CPP = g++ $(CPPOPT)
#CPP = icc $(CPPOPT)
OPTLIB = -lgf2x

LINKOPT += -L../lib -lMTToolBox -lUnitTest++ -lntl $(OPTLIB)
COMMON_OBJS = tinymt32.o tinymt64.o
OBJS1 = ${COMMON_OBJS} test_period.o test_search.o test_equidist.o
OBJS2 = ${COMMON_OBJS} test_temper.o test_search_temper.o

DEPENDS = test_generator.hpp tinymt32.h tinymt64.h

main2: main.cpp ${OBJS2} ${DEPENDS}
	${CPP} -o $@ main.cpp ${OBJS2} ${LINKOPT}

main1: main.cpp ${OBJS1} ${DEPENDS}
	${CPP} -o $@ main.cpp ${OBJS1} ${LINKOPT}

test_rtiny: test_rtiny.cpp test_generator.hpp
	${CPP} -o $@ $< ${LINKOPT}

test_rlittle: test_rlittle.cpp test_generator.hpp
	${CPP} -o $@ $< ${LINKOPT}

doc:
	doxygen doxygen.cfg

.cpp.o:
	${CPP} -c $<

.c.o:
	${CC} -c $<

clean:
	rm -rf *.o *~ *.dSYM
