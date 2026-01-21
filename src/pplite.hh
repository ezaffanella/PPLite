/* PPLite: a lightweight library for convex polyhedra derived from PPL.
   Copyright (C) 2018-2026 Enea Zaffanella <enea.zaffanella@unipr.it>

This file is part of PPLite.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "pplite-config.h"
#include "globals.hh"
#include "clock.hh"
#include "Local_Stats.hh"
#include "Low_Level_Stats.hh"
#include "utils.hh"
#include "support_utils.hh"
#include "ascii_dump_load.hh"
#include "ascii_dump.hh"
#include "Bits.hh"
#include "Index_Partition.hh"
#include "Integer_fwd.hh"
#include "Integer.hh"
#include "Rational_fwd.hh"
#include "Rational.hh"
#include "Var.hh"
#include "Output_Function.hh"
#include "Linear_Expr.hh"
#include "Affine_Expr.hh"
#include "Itv.hh"
#include "BBox.hh"
#include "Con.hh"
#include "Gen.hh"
#include "mater_iterator.hh"
#include "Scalar_Prod.hh"
#include "Sat.hh"
#include "Poly_Rel.hh"
#include "Poly.hh"
#include "Poly_templ.hh"
#include "Poly_min.hh"
#include "Poly_widen.hh"
#include "Poly_Stats.hh"
#include "B_Poly.hh"
#include "F_Poly.hh"
#include "Two_Poly.hh"
#include "U_Poly.hh"
#include "PolySet.hh"
#include "PolySet_templ.hh"
#include "Abs_Poly.hh"
#include "Abs_Poly_Adapter.hh"
#include "Dyn_Poly.hh"
#include "memory_in_bytes.hh"
