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
#include "Poly_Rel.hh"
#include <iostream>
#include <cassert>

namespace pplite {

void
Poly_Con_Rel::ascii_dump(std::ostream& s) const {
  auto f = flags;
  if (f == NOTHING) {
    s << "NOTHING";
    return;
  }

  while (true) {
    if (implies(f, IS_DISJOINT)) {
      s << "IS_DISJOINT";
      f &= ~IS_DISJOINT;
    }
    else if (implies(f, STRICTLY_INTERSECTS)) {
      s << "STRICTLY_INTERSECTS";
      f &= ~STRICTLY_INTERSECTS;
    }
    else if (implies(f, IS_INCLUDED)) {
      s << "IS_INCLUDED";
      f &= ~IS_INCLUDED;
    }
    else if (implies(f, SATURATES)) {
      s << "SATURATES";
      f &= ~SATURATES;
    }
    if (f != NOTHING) {
      s << " & ";
    }
    else {
      break;
    }
  }
}

void
Poly_Gen_Rel::ascii_dump(std::ostream& s) const {
  if (impl() == NOTHING)
    s << "NOTHING";
  else {
    assert(impl() == SUBSUMES);
    s << "SUBSUMES";
  }
}

} // namespace pplite
