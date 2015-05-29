OBJS = tinymt32.o tinymt64.o
OBJS1 = test_period.o test_search.o test_equidist.o
OBJS2 = test_temper.o test_search_temper.o \
        test_linearity.o test_mt.o test_primefactors.o
CXXFLAGS = $(CXXFLAGS) -I../include

test: test1 test2

test1:main.cpp $(OBJS) $(OBJS1)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(OBJS1)

test2:main.cpp $(OBJS) $(OBJS2)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(OBJS2)