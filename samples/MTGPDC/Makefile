#for GNU make

#DDEBUG = -O0 -g -ggdb -DDEBUG=1

# uncomment next line if you installed NTL with gf2x or gmp
OPTLIB = -lgf2x
#OPTLIB = -lgmp

INC = -I../../include

CC = gcc -Wall -O3 -std=c99 -Wmissing-prototypes -Wextra $(DDEBUG)

CPPOPT = -Wall -Wextra -O3 $(M64) $(INC) \
	$(DDEBUG) -D__STDC_CONSTANT_MACROS -D__STDC_FORMAT_MACROS
CPP = g++ $(CPPOPT)
#CPP = icc $(CPPOPT)

LINKOPT += -L../../lib -lntl -lMTToolBox -lcrypto $(OPTLIB)

MTGP32dc:  MTGP32dc.cpp MTGP32search.hpp mtgp_param.h \
	parse_opt.hpp parse_opt.o
	${CPP} -DUSE_SHA -o $@ MTGP32dc.cpp parse_opt.o ${LINKOPT}

test: test.cpp test_linearity.cpp
	${CPP} -DUSE_SHA -o $@ test.cpp test_linearity.cpp -lUnitTest++

.cpp.o:
	${CPP} -c $<

.c.o:
	${CC} -c $<

clean:
	rm -rf *.o *~ *.dSYM
