/* PPLite: a lightweight library for convex polyhedra derived from PPL.
   Copyright (C) 2001-2010 Roberto <bagnara@cs.unipr.it>
   Copyright (C) 2010-2018 BUGSENG srl (http://bugseng.com)
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
test02() {
  Var A(0); Var B(1);
  Var E(4);

  Poly ph(7, Topol::NNC);
  ph.add_con(E < 0);
  ph.add_con(B <= 1);
  ph.add_con(B >= 0);
  ph.add_con(A >= 20);
  ph.minimize();

  Dims pfunc(7, not_a_dim());
  pfunc[2] = 0;
  pfunc[6] = 1;

  ph.map_space_dims(pfunc);
  ph.minimize();

  return (ph.check_inv());
}

bool
test03() {
  Var X(29);
  Poly ph(39, Topol::NNC);
  ph.add_con(Var(25) == 1);
  ph.add_con(Var(26) == -1);
  ph.add_con(X <= 1);
  ph.add_con(Var(23) <= 20);
  ph.add_con(X >= 0);
  ph.minimize();

  ph.add_con(X == 1);

  ph.minimize();

  return (ph.check_inv());
}

bool
test04() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Topol::NNC);
  ph1.add_con(A >= 0);
  ph1.add_con(B > 0);
  ph1.add_con(A <= 2);
  ph1.add_gen(ray(A));
  // Here 1 >= 0 comes back in sk.
  ph1.minimize();
  ph1.ascii_dump(nout);
  // Pos is non-strict and should be promoted in 1 > 0.
  ph1.topological_closure_assign();

  Poly ph2(2, Topol::NNC);
  ph2.add_con(A >= 0);
  ph2.add_con(B >= 0);

  ph1.ascii_dump(nout);
  bool ok = (ph1 == ph2);
  return ok;
}

BEGIN_MAIN
/* test01 moved into Poly_Impl */
  DO_TEST(test02);
  DO_TEST(test03);
  DO_TEST(test04);
END_MAIN
