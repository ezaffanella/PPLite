dnl PPLite: a lightweight library for convex polyhedra derived from PPL.
dnl Copyright (C) 2018-2023 Enea Zaffanella <enea.zaffanella@unipr.it>
dnl
dnl This file is part of PPLite.
dnl
dnl This program is free software: you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation, either version 3 of the License, or
dnl (at your option) any later version.
dnl
dnl This program is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.
dnl
dnl You should have received a copy of the GNU General Public License
dnl along with this program.  If not, see <http://www.gnu.org/licenses/>.

dnl A function to check for the existence and usability of mpfr.

AC_DEFUN([AC_CHECK_MPFR],
[
AC_ARG_WITH(mpfr,
  AS_HELP_STRING([--with-mpfr=DIR],
		 [search for libmpfr in DIR/include and DIR/lib]))

AC_ARG_WITH(mpfr-include,
  AS_HELP_STRING([--with-mpfr-include=DIR],
		 [search for libmpfr headers in DIR]))

AC_ARG_WITH(mpfr-lib,
  AS_HELP_STRING([--with-mpfr-lib=DIR],
		 [search for libmpfr library objects in DIR]))

if test -n "$with_mpfr"
then
  mpfr_include_options="-I$with_mpfr/include"
  mpfr_library_paths="$with_mpfr/lib"
  mpfr_library_options="-L$mpfr_library_paths"
fi

if test -n "$with_mpfr_include"
then
  mpfr_include_options="-I$with_mpfr_include"
fi

if test -n "$with_mpfr_lib"
then
  mpfr_library_paths="$with_mpfr_lib"
  mpfr_library_options="-L$mpfr_library_paths"
fi

mpfr_libs="-lmpfr"

mpfr_library_options="$mpfr_library_options $mpfr_libs"

ac_save_CPPFLAGS="$CPPFLAGS"
CPPFLAGS="$CPPFLAGS $mpfr_include_options"
ac_save_LIBS="$LIBS"
LIBS="$LIBS $mpfr_library_options"

AC_LANG_PUSH(C)

AC_MSG_CHECKING([for the MPFR library])
AC_RUN_IFELSE([AC_LANG_SOURCE([[
#include <mpfr.h>
int main() {
  mpfr_t x;
  mpfr_init2(x, 200);
  mpfr_clear(x);
  return 0;
}
]])],
  AC_MSG_RESULT(yes)
  ac_cv_have_mpfr=yes,
  AC_MSG_RESULT(no)
  ac_cv_have_mpfr=no)

AC_LANG_POP(C)
LIBS="$ac_save_LIBS"
CPPFLAGS="$ac_save_CPPFLAGS"

])
