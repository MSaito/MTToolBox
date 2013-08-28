#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([MTToolBox], [0.2], [saito@math.sci.hiroshima-u.ac.jp])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([test/main.cpp])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP

# Checks for libraries.
AC_PROG_RANLIB
# FIXME: Replace `main' with a function in `-lcrypto':
AC_CHECK_LIB([crypto], [SHA1_Update])
# FIXME: Replace `main' with a function in `-lgf2x':
AC_CHECK_LIB([gf2x], [gf2x_mul_1_n])
AC_CHECK_LIB([gmp], [__gmpz_init])

# Checks for header files.
AC_CHECK_HEADERS([inttypes.h stdint.h stdlib.h unistd.h, NTL/version.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T

# Checks for library functions.
AC_CHECK_FUNCS([memset strtol strtoull])

AC_CONFIG_FILES([lib/Makefile
                 samples/MTGPDC/Makefile
                 samples/TinyMTDC/Makefile
                 samples/sfmtdc/Makefile
                 test/Makefile])
AC_OUTPUT(Makefile)
