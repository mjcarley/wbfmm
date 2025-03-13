AC_INIT([wbfmm],[0.1])
AC_CONFIG_SRCDIR([src/wbfmm.h])
AC_CONFIG_MACRO_DIR([m4])

## this optimization is pretty much essential for WBFMM to be useful
: ${CFLAGS=-O3 -g} 

WBFMM_MAJOR_VERSION=0
WBFMM_MINOR_VERSION=1
WBFMM_MICRO_VERSION=1
WBFMM_INTERFACE_AGE=0
WBFMM_BINARY_AGE=0
WBFMM_VERSION=$WBFMM_MAJOR_VERSION.$WBFMM_MINOR_VERSION.$WBFMM_MICRO_VERSION

AC_SUBST(WBFMM_MAJOR_VERSION)
AC_SUBST(WBFMM_MINOR_VERSION)
AC_SUBST(WBFMM_MICRO_VERSION)
AC_SUBST(WBFMM_VERSION)

# libtool versioning
LT_RELEASE=$WBFMM_MAJOR_VERSION.$WBFMM_MINOR_VERSION
LT_CURRENT=`expr $WBFMM_MICRO_VERSION - $WBFMM_INTERFACE_AGE`
LT_REVISION=$WBFMM_INTERFACE_AGE
LT_AGE=`expr $WBFMM_BINARY_AGE - $WBFMM_INTERFACE_AGE`
AC_SUBST(LT_RELEASE)
AC_SUBST(LT_CURRENT)
AC_SUBST(LT_REVISION)
AC_SUBST(LT_AGE)

# For automake.
VERSION=$WBFMM_VERSION
PACKAGE=wbfmm

AC_SUBST(PACKAGE)
AC_SUBST(VERSION)

AM_INIT_AUTOMAKE(subdir-objects)

##$PACKAGE, $VERSION, no-define)

# Specify a configuration file
AM_CONFIG_HEADER(config.h)

AC_DEFINE_UNQUOTED(WBFMM_MAJOR_VERSION, $WBFMM_MAJOR_VERSION, [Major version])
AC_DEFINE_UNQUOTED(WBFMM_MINOR_VERSION, $WBFMM_MINOR_VERSION, [Minor version])
AC_DEFINE_UNQUOTED(WBFMM_MICRO_VERSION, $WBFMM_MICRO_VERSION, [Micro version])
AC_DEFINE_UNQUOTED(WBFMM_INTERFACE_AGE, $WBFMM_INTERFACE_AGE, [Interface age])
AC_DEFINE_UNQUOTED(WBFMM_BINARY_AGE, $WBFMM_BINARY_AGE, [Binary age])

dnl Initialize libtool
LT_INIT

dnl Initialize maintainer mode
#AM_MAINTAINER_MODE

AC_PROG_CC
AM_PROG_CC_C_O

if test x$GCC = xyes ; then
  CFLAGS="$CFLAGS -Wall -Werror-implicit-function-declaration -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations"
fi

AC_PROG_AWK
AC_SEARCH_LIBS([strerror],[cposix])
AC_HEADER_STDC

PKG_CHECK_MODULES([GLIB], [glib-2.0], [],
[
  echo \#################################################################
  echo \# 
  echo \# GLIB-2.0 not installed. You can install it using your 
  echo \# distribution's package manager or from:
  echo \#
  echo \# https://gitlab.gnome.org/GNOME/glib/
  echo \#
  echo \#################################################################
  exit
])

PKG_CHECK_MODULES([BLASWRAP], [blaswrap], [],
[
  echo \#################################################################
  echo \# 
  echo \# BLASWRAP not installed. You can install it using your 
  echo \# distribution's package manager or from:
  echo \#
  echo \# https://github.com/mjcarley/blaswrap
  echo \#
  echo \#################################################################
  exit
])

AX_GCC_X86_CPU_SUPPORTS(avx, hasavx=yes) 
AX_GCC_X86_CPU_SUPPORTS(avx2, hasavx2=yes) 
AX_GCC_X86_CPU_SUPPORTS(fma, hasfma=yes) 

if test x$hasavx = xyes; then
   SIMD_FLAGS="$SIMD_FLAGS -mavx"
fi
if test x$hasfma = xyes; then
   SIMD_FLAGS="$SIMD_FLAGS -mfma"
fi

AC_CHECK_LIB(m, cos)
AC_OPENMP
AC_PROG_F77
AC_F77_LIBRARY_LDFLAGS
AC_F77_DUMMY_MAIN
AC_F77_WRAPPERS
##ACX_BLAS

#DX_INIT_DOXYGEN(wbfmm, doxy/wbfmm.dxy, doc)
#DX_DOXYGEN_FEATURE(ON)
#DX_DOT_FEATURE(ON)
#DX_HTML_FEATURE(ON)
#DX_CHM_FEATURE(OFF)
#DX_CHI_FEATURE(OFF)
#DX_MAN_FEATURE(OFF)
#DX_RTF_FEATURE(OFF)
#DX_XML_FEATURE(OFF)
#DX_PDF_FEATURE(OFF)
#DX_PS_FEATURE(OFF)

CFLAGS="$CFLAGS $OPENMP_CFLAGS $SIMD_FLAGS $GLIB_CFLAGS $BLASWRAP_CFLAGS"
LIBS="$GLIB_LIBS $BLASWRAP_LIBS $LAPACK_LIBS $LIBS $FLIBS"

AC_SUBST(CFLAGS)
AC_SUBST(LIBS)
AC_SUBST(LDFLAGS)
#AC_SUBST(DX_RULES)

AC_CONFIG_FILES([
Makefile
wbfmm.pc
src/Makefile
tools/Makefile
])

AC_OUTPUT

