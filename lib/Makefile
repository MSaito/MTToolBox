WARN = -Wall -Wextra
OPTI = -O3
CXX = g++
AR = ar
INCLUDE = -I../include/
OBJS = period.o AlgorithmPrimitivity.o
CXXFLAGS = -D__STDC_CONSTANT_MACROS -D__STDC_FORMAT_MACROS $(OPTI) \
$(WARN) $(STD) $(INCLUDE)

libMTToolBox.a: $(OBJS)
	$(AR) rv $@ $(OBJS)

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $<