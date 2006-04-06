#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(mhdcl, 0.1b, mbaptist@fc.up.pt)

AC_CONFIG_SRCDIR([lss.C])

AC_PROG_INSTALL
AC_PROG_RANLIB

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP

# Checks for libraries.
AC_CHECK_LIB([fftw3],[fftw_destroy_plan])

# Checks for header files.


# Check for cat++
AC_ARG_WITH(cat,
    [ --with-cat=<path> prefix of MySQL installation. e.g. /usr/local or /usr],
    [CAT_PREFIX=$with_cat],
    AC_MSG_ERROR([You must call configure with the --with-cat option.
    This tells configure where to find the cat++ library and headers.
    e.g. --with-cat=/usr/local or --with-cat=/usr])
)
AC_SUBST(CAT_PREFIX)
CAT_LIBS="-L${CAT_PREFIX}"
CAT_CFLAGS="-I${CAT_PREFIX}"
AC_SUBST(CAT_LIBS)
AC_SUBST(CAT_CFLAGS)

# Check for lass
AC_ARG_WITH(lass,
    [ --with-lass=<path> prefix of lass installation. e.g. /usr/local or /usr],
    [LASS_PREFIX=$with_lass],
    AC_MSG_ERROR([You must call configure with the --with-lass option.
    This tells configure where to find the lass library and headers.
    e.g. --with-lass=/usr/local or --with-lass=/usr])
)
AC_SUBST(LASS_PREFIX)
LASS_CFLAGS="-I${LASS_PREFIX}"
AC_SUBST(LASS_CFLAGS)

# Check for goops
AC_ARG_WITH(goops,
    [ --with-goops=<path> prefix of cgsolver installation. e.g. /usr/local or /usr],
    [GOOPS_PREFIX=$with_goops],
    AC_MSG_ERROR([You must call configure with the --with-goops option.
    This tells configure where to find the goops library and headers.
    e.g. --with-goops=/usr/local or --with-goops=/usr])
)
AC_SUBST(GOOPS_PREFIX)
GOOPS_LIBS="-L${GOOPS_PREFIX}"
GOOPS_CFLAGS="-I${GOOPS_PREFIX}"
AC_SUBST(GOOPS_LIBS)
AC_SUBST(GOOPS_CFLAGS)



# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE

# Checks for library functions.

AC_CONFIG_FILES([Makefile vzdeigen/Makefile])
AC_OUTPUT

