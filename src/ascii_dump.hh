/* PPLite: a lightweight library for convex polyhedra derived from PPL.
   Copyright (C) 2018-2024 Enea Zaffanella <enea.zaffanella@unipr.it>

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

#ifndef pplite_ascii_dump_hh
#define pplite_ascii_dump_hh 1

#include "Integer_fwd.hh"
#include "Rational_fwd.hh"

// No way to forward declare vectors (Cons, Gens, NS_Rows, ...).
#include "Bits.hh"
#include "Con.hh"
#include "Gen.hh"

namespace pplite {

// Forward declarations.
class Linear_Expr;
class Affine_Expr;
class Sat;
template <bool> class Box;
class Poly;
struct Poly_Impl;
class Poly_Con_Rel;
class Poly_Gen_Rel;

// Note: the following ascii_dump() functions are provided for debugging,
// mainly to simplify visual inspection of pplite class objects.
// They are *meant* to be pure declarations (defined in ascii_dump.cc);
// in particular, they are meant to be non-template and non-inline,
// since this allows to call them from gdb more easily.

void ascii_dump(const Bits&);
void ascii_dump(const NS_Rows&);
void ascii_dump(const Integer&);
void ascii_dump(const Rational&);
void ascii_dump(const Linear_Expr&);
void ascii_dump(const Affine_Expr&);
void ascii_dump(const Con&);
void ascii_dump(const Cons&);
void ascii_dump(const Gen&);
void ascii_dump(const Gens&);
void ascii_dump(const Sat&);
template <bool KVI>
void ascii_dump(const Box<KVI>&);
void ascii_dump(const Poly&);
void ascii_dump(const Poly_Impl&);
void ascii_dump(const Poly_Con_Rel&);
void ascii_dump(const Poly_Gen_Rel&);

} // namespace pplite

#endif // !defined(pplite_ascii_dump_hh)
