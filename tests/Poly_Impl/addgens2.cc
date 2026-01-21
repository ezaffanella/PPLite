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

#include "pplite_test.hh"

bool
test03() {
  Var A(0);
  Var B(1);

  Poly ph(2, Topol::NNC, Spec_Elem::EMPTY);
  ph.add_gen(point(A));
  ph.add_gen(point(B));
  ph.add_gen(closure_point(A + 2*B));
  ph.add_gen(ray(A));
  ph.minimize();

  auto ns_dim = ph.impl().cs.ns_rows.size();
  bool ok = (ns_dim == 0);

  ph.add_gen(point(2*A + 2*B));
  ph.minimize();
  ns_dim = ph.impl().cs.ns_rows.size();
  ok = ok && (ns_dim == 2);

  ph.add_gen(point(2*B));
  ph.minimize();
  ns_dim = ph.impl().cs.ns_rows.size();
  ok = ok && (ns_dim == 1);

  ph.add_gen(closure_point(3*B));
   ph.minimize();
  ns_dim = ph.impl().cs.ns_rows.size();
  ok = ok && (ns_dim == 0);

  return ok;
}

BEGIN_MAIN
  DO_TEST(test03);
END_MAIN
