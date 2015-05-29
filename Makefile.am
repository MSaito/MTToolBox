EXTRA_DIST = NEWS.jp.utf8.txt README.jp.utf8.txt include test doxygen slide

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = lib samples include

AM_CFLAGS = -O3 -Wall -Wextra -Wno-unused-parameter -std=c99
AM_LIBS = -lntl

doc:
	$(MAKE) -C doxygen doc