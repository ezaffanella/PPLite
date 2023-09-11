# Public macros for the TeX Live (TL) tree.
# Copyright (C) 1995-2009 Karl Berry <tex-live@tug.org>
# Copyright (C) 2009-2012 Peter Breitenlohner <tex-live@tug.org>
#
# This file is free software; the copyright holders
# give unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.


# Use libtool for linking.  This allows testing properties of libraries
# that can be either (1) uninstalled libtool libraries already built
# when this configure runs, or (2) installed libraries -- libtool or not.

AC_DEFUN([AC_USE_LIBTOOL],
[## $0: Generate a libtool script for use in configure tests
AC_PROVIDE_IFELSE([LT_INIT], ,
                  [m4_fatal([$0: requires libtool])])[]dnl
LT_OUTPUT
m4_append([AC_LANG(C)],
[ac_link="./libtool --mode=link --tag=CC $ac_link"
])[]dnl
AC_PROVIDE_IFELSE([AC_PROG_CXX],
[m4_append([AC_LANG(C++)],
[ac_link="./libtool --mode=link --tag=CXX $ac_link"
])])[]dnl
AC_LANG(_AC_LANG)[]dnl
]) # AC_USE_LIBTOOL
