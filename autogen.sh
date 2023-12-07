#! /bin/sh

mkdir -p m4

aclocal -I . \
&& automake --add-missing \
&& autoconf
