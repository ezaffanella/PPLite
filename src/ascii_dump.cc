/* PPLite: a lightweight library for convex polyhedra derived from PPL.
   Copyright (C) 2001-2010 Roberto Bagnara <bagnara@cs.unipr.it>
   Copyright (C) 2010-2018 BUGSENG srl (http://bugseng.com)
   Copyright (C) 2018-2023 Enea Zaffanella <enea.zaffanella@unipr.it>

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

#include "ascii_dump.hh"
#include "Bits.hh"
#include "Integer.hh"
#include "Rational.hh"
#include "Linear_Expr.hh"
#include "Affine_Expr.hh"
#include "Con.hh"
#include "Gen.hh"
#include "Sat.hh"
#include "BBox.hh"
#include "Poly.hh"
#include "Poly_templ.hh"
#include "Poly_Rel.hh"

#include <iostream>
#include <vector>

namespace pplite {

namespace detail {

template <typename T>
void ascii_dump_impl(const T& t) {
  t.ascii_dump(std::cerr);
}

template <typename T>
void ascii_dump_impl(const std::vector<T>& vt) {
  for (const auto& t : vt)
    ascii_dump_impl(t);
}

} // namespace detail

// Note: the following ascii_dump() functions are provided for debugging,
// mainly to simplify visual inspection of pplite class objects.
// They are *meant* to be non-template and non-inline,
// since this allows to call them from gdb more easily.

using detail::ascii_dump_impl;

void ascii_dump(const Bits& x) { ascii_dump_impl(x); }
void ascii_dump(const NS_Rows& x) { ascii_dump_impl(x); }
void ascii_dump(const Integer& x) { ascii_dump_impl(x); }
void ascii_dump(const Rational& x) { ascii_dump_impl(x); }
void ascii_dump(const Linear_Expr& x) { ascii_dump_impl(x); }
void ascii_dump(const Affine_Expr& x) { ascii_dump_impl(x); }
void ascii_dump(const Con& x) { ascii_dump_impl(x); }
void ascii_dump(const Cons& x) { ascii_dump_impl(x); }
void ascii_dump(const Gen& x) { ascii_dump_impl(x); }
void ascii_dump(const Gens& x) { ascii_dump_impl(x); }
void ascii_dump(const Sat& x) { ascii_dump_impl(x); }
template <bool KVI>
void ascii_dump(const Box<KVI>& x) { ascii_dump_impl(x); }
void ascii_dump(const Poly& x) { ascii_dump_impl(x); }
void ascii_dump(const Poly_Impl& x) { ascii_dump_impl(x); }
void ascii_dump(const Poly_Impl::Sys<Cons>& x) { ascii_dump_impl(x); }
void ascii_dump(const Poly_Impl::Sys<Gens>& x) { ascii_dump_impl(x); }
void ascii_dump(const Poly_Con_Rel& x) { ascii_dump_impl(x); }
void ascii_dump(const Poly_Gen_Rel& x) { ascii_dump_impl(x); }

} // namespace pplite

