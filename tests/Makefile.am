## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/src -I$(includedir)

DEPS = $(top_builddir)/src/libwbfmm.la
LDADD =	$(top_builddir)/src/.libs/libwbfmm.a

EXTRA_DIST = \
	wbfmm-tests.c.in wbfmm-tree-calc.c.in \
	wbfmm-direct-calc.c.in \
	wbfmm-laplace-tests.c.in wbfmm-laplace-tree-calc.c.in \
	make-source make-source-float

bin_PROGRAMS = \
	wbfmm-tests \
	wbfmm-tests-avx \
	wbfmm-tree-calc \
	wbfmm-tree-calc-avx \
	wbfmm-tests-float \
	wbfmm-tree-calc-float \
	wbfmm-direct-calc \
	wbfmm-direct-calc-float \
	wbfmm-laplace-tests \
	wbfmm-laplace-tests-avx \
	wbfmm-laplace-tests-float \
	wbfmm-laplace-tree-calc \
	wbfmm-laplace-tree-calc-avx \
	wbfmm-laplace-tree-calc-float \
	wbfmm-laplace-direct-calc \
	wbfmm-laplace-direct-calc-float

BUILT_SOURCES = \
	wbfmm-tests.c \
	wbfmm-tree-calc.c \
	wbfmm-direct-calc.c \
	wbfmm-tests-float.c \
	wbfmm-tree-calc-float.c \
	wbfmm-direct-calc-float.c \
	wbfmm-laplace-tests.c \
	wbfmm-laplace-tests-float.c \
	wbfmm-laplace-tree-calc.c \
	wbfmm-laplace-tree-calc-float.c \
	wbfmm-laplace-direct-calc.c \
	wbfmm-laplace-direct-calc-float.c

CLEANFILES = \
	wbfmm-tests.c wbfmm-tests-float.c \
	wbfmm-tree-calc.c wbfmm-tree-calc-float.c \
	wbfmm-direct-calc.c wbfmm-direct-calc-float.c \
	wbfmm-laplace-tests.c wbfmm-laplace-tests-float.c \
	wbfmm-laplace-tree-calc.c wbfmm-laplace-tree-calc-float.c \
	wbfmm-laplace-direct-calc.c wbfmm-laplace-direct-calc-float.c 

wbfmm-tests.c: wbfmm-tests.c.in
	./make-source wbfmm-tests

wbfmm-tests-float.c: wbfmm-tests.c.in
	./make-source-float wbfmm-tests

wbfmm-tree-calc.c: wbfmm-tree-calc.c.in
	./make-source wbfmm-tree-calc

wbfmm-tree-calc-float.c: wbfmm-tree-calc.c.in
	./make-source-float wbfmm-tree-calc

wbfmm-direct-calc.c: wbfmm-direct-calc.c.in
	./make-source wbfmm-direct-calc

wbfmm-direct-calc-float.c: wbfmm-direct-calc.c.in
	./make-source-float wbfmm-direct-calc

wbfmm-laplace-tests.c: wbfmm-laplace-tests.c.in
	./make-source wbfmm-laplace-tests

wbfmm-laplace-tests-float.c: wbfmm-laplace-tests.c.in
	./make-source-float wbfmm-laplace-tests

wbfmm-laplace-tree-calc.c: wbfmm-laplace-tree-calc.c.in
	./make-source wbfmm-laplace-tree-calc

wbfmm-laplace-tree-calc-float.c: wbfmm-laplace-tree-calc.c.in
	./make-source-float wbfmm-laplace-tree-calc

wbfmm-laplace-direct-calc.c: wbfmm-laplace-direct-calc.c.in
	./make-source wbfmm-laplace-direct-calc

wbfmm-laplace-direct-calc-float.c: wbfmm-laplace-direct-calc.c.in
	./make-source-float wbfmm-laplace-direct-calc

wbfmm_tests_LDADD = \
	$(top_builddir)/src/.libs/libwbfmm.a

wbfmm_tests_avx_SOURCES = wbfmm-tests.c
wbfmm_tests_avx_CFLAGS = -DWBFMM_USE_AVX $(AM_CFLAGS)

wbfmm_tree_calc_SOURCES = wbfmm-tree-calc.c

wbfmm_tree_calc_avx_SOURCES = wbfmm-tree-calc.c
wbfmm_tree_calc_avx_CFLAGS = -DWBFMM_USE_AVX $(AM_CFLAGS)

wbfmm_laplace_tests_avx_SOURCES = wbfmm-laplace-tests.c
wbfmm_laplace_tests_avx_CFLAGS = -DWBFMM_USE_AVX $(AM_CFLAGS)

wbfmm_laplace_tree_calc_avx_SOURCES = wbfmm-laplace-tree-calc.c

wbfmm_laplace_tree_calc_avx_CFLAGS = -DWBFMM_USE_AVX $(AM_CFLAGS)
