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

#include "pplite_test.hh"

bool
test11() {
  Var A(0);
  Var B(1);

  Poly ph(2, Spec_Elem::UNIVERSE, Topol::NNC);
  ph.add_con(B <= 2);
  ph.add_con(B > 0);
  ph.add_con(A + B >= 1);
  ph.add_con(A - B >= -1);

  ph.add_gen(point(2*A));

  Poly known_result(2, Spec_Elem::EMPTY, Topol::NNC);
  known_result.add_gen(point(B));
  known_result.add_gen(point(A + 2*B));
  known_result.add_gen(closure_point(A));
  known_result.add_gen(ray(A));
  known_result.add_gen(point(2*A));

  bool ok = (ph == known_result)
    && ph.impl().cs.ns_rows.size() == 2;
  return ok;
}

bool
test13() {
  Var A(0);
  Var B(1);
  Var C(2);
  Var D(3);
  Var E(4);
  Var F(5);

  Poly ph1(6, Topol::NNC);
  ph1.add_con(A == 1);
  ph1.add_con(E -10*F == 0);
  ph1.add_con(E < 0);
  ph1.add_con(B <= 1);
  ph1.add_con(B >= 0);
  Index_Set ns;
  ns.set(1);
  ns.set(2);
  ph1.impl().cs_pending.ns_rows.push_back(ns);

  ph1.ascii_dump(nout);
  ph1.minimize();

  bool ok = (ph1.check_inv());
  return ok;
}

BEGIN_MAIN
  DO_TEST(test11);
  DO_TEST(test13);
END_MAIN

