# PPLite: a lightweight library for convex polyhedra derived from PPL.
# Copyright (C) 2018-2024 Enea Zaffanella <enea.zaffanella@unipr.it>
#
# This file is part of PPLite.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

AM_CPPFLAGS = -I$(top_builddir) -I$(top_srcdir)/src @extra_includes@

if USE_FLINT_INTEGERS

INTEGER_HH_FILES = FLINT_Integer.hh
INTEGER_CC_FILES = FLINT_Integer.cc
RATIONAL_HH_FILES = FLINT_Rational.hh
RATIONAL_CC_FILES = FLINT_Rational.cc

else USE_FLINT_INTEGERS_FALSE

INTEGER_HH_FILES = GMP_Integer.hh
INTEGER_CC_FILES = GMP_Integer.cc
RATIONAL_HH_FILES = GMP_Rational.hh
RATIONAL_CC_FILES = GMP_Rational.cc

endif USE_FLINT_INTEGERS_FALSE

PPLITE_HH_FILES = \
pplite.hh \
ascii_dump_load.hh \
ascii_dump.hh \
clock.hh \
globals.hh \
mater_iterator.hh \
memory_in_bytes.hh \
support_utils.hh \
utils.hh \
Abs_Poly.hh \
Abs_Poly_Adapter.hh \
Affine_Expr.hh \
BBox.hh \
BBox_impl.hh \
Bits.hh \
B_Poly.hh \
Con.hh \
Dyn_Poly.hh \
F_Poly.hh \
Gen.hh \
Index_Partition.hh \
Integer_fwd.hh \
Integer.hh \
$(INTEGER_HH_FILES) \
Itv.hh \
Linear_Expr.hh \
Local_Stats.hh \
Low_Level_Stats.hh \
Output_Function.hh \
Poly.hh \
Poly_min.hh \
Poly_templ.hh \
Poly_widen.hh \
Poly_Rel.hh \
Poly_Stats.hh \
Rational_fwd.hh \
Rational.hh \
$(RATIONAL_HH_FILES) \
Sat.hh \
Scalar_Prod.hh \
Two_Poly.hh \
U_Poly.hh \
PolySet.hh \
PolySet_templ.hh \
Var.hh

PPLITE_CC_FILES = \
ascii_dump.cc \
Abs_Poly_Adapter.cc \
Affine_Expr.cc \
Bits.cc \
B_Poly.cc \
Con.cc \
F_Poly.cc \
Gen.cc \
globals.cc \
$(INTEGER_CC_FILES) \
Linear_Expr.cc \
Low_Level_Stats.cc \
Poly.cc \
Poly_min.cc \
Poly_conv.cc \
Poly_simpl.cc \
Poly_split.cc \
Poly_widen.cc \
Poly_Rel.cc \
Poly_Stats.cc \
$(RATIONAL_CC_FILES) \
Sat.cc \
U_Poly.cc \
PolySet.cc \
Var.cc

lib_LTLIBRARIES = libpplite.la

libpplite_la_SOURCES = \
$(PPLITE_HH_FILES) \
$(PPLITE_CC_FILES)

libpplite_la_LIBADD = @extra_libraries@

# PPLite version 0.7 --> libtool version info 1:0:0
# PPLite version 0.8 --> libtool version info 2:0:0
# PPLite version 0.9 --> libtool version info 3:0:0
# PPLite version 0.10 --> libtool version info 4:0:0
# PPLite version 0.12 --> libtool version info 5:0:0
libpplite_la_LDFLAGS = -no-undefined -version-info 5:0:0

pkginclude_HEADERS = \
$(top_builddir)/pplite-config.h \
$(PPLITE_HH_FILES)
