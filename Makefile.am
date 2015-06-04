EXTRA_DIST = NEWS.jp.utf8.txt README.jp.utf8.txt include test doxygen slide

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = lib samples include test

AM_CFLAGS = -O2 -Wall -Wextra -std=c99
AM_CXXFLAGS = -O2 -Wall -Wextra -Wno-unused-parameter -std=c++11
AM_LIBS = -lntl

doc:
	$(MAKE) -C doxygen doc