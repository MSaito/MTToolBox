#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([MTToolBox], [0.2.1], [saito@manieth.com])
AC_PACKAGE_NAME('MTToolBox')
AC_PACKAGE_TARNAME('MTToolBox')
AC_CONFIG_SRCDIR([lib/period.cpp])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIRS([m4])
AM_INIT_AUTOMAKE
#LT_INIT
AC_PROG_RANLIB
AC_SUBST(PACKAGE_NAME, 'MTToolBox')

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([inttypes.h limits.h stdint.h stdlib.h string.h unistd.h])
#AC_CHECK_

# Checks for typedefs, structures, and compiler characteristics.
#AC_CHECK_HEADER_STDBOOL
#AC_C_INLINE
#AC_TYPE_INT16_T
#AC_TYPE_INT32_T
#AC_TYPE_INT64_T
#AC_TYPE_UINT16_T
#AC_TYPE_UINT32_T
#AC_TYPE_UINT64_T

AC_CHECK_LIB(gf2x, gf2x_mul)
AC_CHECK_LIB(gmp, __gmpz_init)

# Checks for library functions.
#AC_CHECK_FUNCS([memset strtol strtoull])

AC_CONFIG_FILES([Makefile include/Makefile lib/Makefile samples/Makefile \
			  samples/MTDC/Makefile samples/MTGPDC/Makefile \
			  samples/XORSHIFT/Makefile samples/RMT/Makefile \
			  samples/TinyMTDC/Makefile test/Makefile])
AC_OUTPUT
