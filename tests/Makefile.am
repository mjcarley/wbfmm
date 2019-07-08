## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/src -I$(includedir)

DEPS = $(top_builddir)/src/libwbfmm.la

EXTRA_DIST = wbfmm-test.c.in wbfmm-tree-calc.c.in

bin_PROGRAMS = \
	wbfmm-tests \
	wbfmm-tree-calc \
	wbfmm-tests-float \
	wbfmm-tree-calc-float \
	wbfmm-direct-calc \
	wbfmm-direct-calc-float

BUILT_SOURCES = \
	wbfmm-tests.c \
	wbfmm-tree-calc.c \
	wbfmm-direct-calc.c \
	wbfmm-tests-float.c \
	wbfmm-tree-calc-float.c
	wbfmm-direct-calc-float.c

CLEANFILES = wbfmm-tests.c wbfmm-tests-float.c \
	wbfmm-tree-calc.c wbfmm-tree-calc-float.c \
	wbfmm-direct-calc.c wbfmm-direct-calc-float.c

wbfmm-tests.c: wbfmm-tests.c.in
	sed -e s/WBFMM_REAL/gdouble/g -e s/_PREC//g -e s/FLOATFMT/lg/g \
	-e s/FLOATEXP/1.16e/g \
	wbfmm-tests.c.in > wbfmm-tests.c

wbfmm-tests-float.c: wbfmm-tests.c.in
	sed -e s/WBFMM_REAL/gfloat/g -e s/_PREC/_f/g -e s/FLOATFMT/g/g \
	-e s/FLOATEXP/1.16e/g \
	wbfmm-tests.c.in > wbfmm-tests-float.c

wbfmm-tree-calc.c: wbfmm-tree-calc.c.in
	sed -e s/WBFMM_REAL/gdouble/g -e s/_PREC//g -e s/FLOATFMT/lg/g \
	-e s/FLOATEXP/1.16e/g \
	wbfmm-tree-calc.c.in > wbfmm-tree-calc.c

wbfmm-tree-calc-float.c: wbfmm-tree-calc.c.in
	sed -e s/WBFMM_REAL/gfloat/g -e s/_PREC/_f/g -e s/FLOATFMT/g/g \
	-e s/FLOATEXP/1.16e/g \
	wbfmm-tree-calc.c.in > wbfmm-tree-calc-float.c

wbfmm-direct-calc.c: wbfmm-direct-calc.c.in
	sed -e s/WBFMM_REAL/gdouble/g -e s/_PREC//g -e s/FLOATFMT/lg/g \
	-e s/FLOATEXP/1.16e/g \
	wbfmm-direct-calc.c.in > wbfmm-direct-calc.c

wbfmm-direct-calc-float.c: wbfmm-direct-calc.c.in
	sed -e s/WBFMM_REAL/gfloat/g -e s/_PREC/_f/g -e s/FLOATFMT/g/g \
	-e s/FLOATEXP/1.16e/g \
	wbfmm-direct-calc.c.in > wbfmm-direct-calc-float.c

wbfmm_tests_LDADD = \
	$(top_builddir)/src/.libs/libwbfmm.a \
	$(top_builddir)/src/.libs/libwbfmm-base.a

wbfmm_tests_float_LDADD = \
	$(top_builddir)/src/.libs/libwbfmm_f.a \
	$(top_builddir)/src/.libs/libwbfmm-base.a

wbfmm_tree_calc_LDADD = \
	$(top_builddir)/src/.libs/libwbfmm.a \
	$(top_builddir)/src/.libs/libwbfmm-base.a

wbfmm_tree_calc_float_LDADD = \
	$(top_builddir)/src/.libs/libwbfmm_f.a \
	$(top_builddir)/src/.libs/libwbfmm-base.a

wbfmm_direct_calc_LDADD = \
	$(top_builddir)/src/.libs/libwbfmm.a \
	$(top_builddir)/src/.libs/libwbfmm-base.a

wbfmm_direct_calc_float_LDADD = \
	$(top_builddir)/src/.libs/libwbfmm_f.a \
	$(top_builddir)/src/.libs/libwbfmm-base.a
