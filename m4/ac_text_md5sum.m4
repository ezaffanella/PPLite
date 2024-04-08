dnl PPLite: a lightweight library for convex polyhedra derived from PPL.
dnl Copyright (C) 2001-2010 Roberto Bagnara <bagnara@cs.unipr.it>
dnl Copyright (C) 2010-2017 BUGSENG srl (http://bugseng.com)
dnl Copyright (C) 2018-2024 Enea Zaffanella <enea.zaffanella@unipr.it>
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

dnl A function to set the command for computing the MD5 checksum of text files.

AC_DEFUN([AC_TEXT_MD5SUM],
[
AC_MSG_CHECKING([for the text md5sum command])
if echo a | (md5sum -t) >/dev/null 2>&1
then
  ac_cv_prog_text_md5sum='md5sum -t'
else
  ac_cv_prog_text_md5sum='md5sum'
fi
AC_MSG_RESULT($ac_cv_prog_text_md5sum)
TEXT_MD5SUM=$ac_cv_prog_text_md5sum
AC_SUBST([TEXT_MD5SUM])
])
