/* PPLite: a lightweight library for convex polyhedra derived from PPL.
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
#include "Affine_Expr.hh"
#include "ascii_dump_load.hh"

#include <iostream>

namespace pplite {

void Affine_Expr::print(std::ostream& os) const {
  if (expr.is_zero()) {
    inhomo.print(os);
    return;
  }
  expr.print(os);
  if (inhomo.is_zero())
    return;
  if (inhomo > 0)
    os << " + ";
  inhomo.print(os);
}

void
Affine_Expr::ascii_dump(std::ostream& s) const {
  s << "expr ";
  expr.ascii_dump(s);
  s << "\n";
  s << "inhomo ";
  inhomo.ascii_dump(s);
}

bool
Affine_Expr::ascii_load(std::istream& s) {
  return ascii_load_string(s, "expr")
    && expr.ascii_load(s)
    && ascii_load_string(s, "inhomo")
    && inhomo.ascii_load(s);
}

} // namespace pplite
