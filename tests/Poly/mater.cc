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
test01() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Topol::NNC);
  ph1.add_con(A >= 0);
  ph1.add_con(A <= 1);
  ph1.add_con(B >= 0);
  ph1.add_con(A + B > 0);
  ph1.add_con(3*A + 2*B > 0);
  ph1.minimize();

  Poly ph2(2, Topol::NNC);
  for (const auto& c : ph1.cons())
    ph2.add_con(c);

  return (ph1 == ph2);
}

bool
test02() {
  Var A(0);
  Var B(1);
  Var C(2);

  Poly ph1(3, Topol::NNC);
  ph1.add_con(A >= 0);
  ph1.add_con(A <= 1);
  ph1.add_con(B >= 0);
  ph1.add_con(A + B > 0);
  ph1.add_con(3*A + 2*B > 0);
  ph1.add_con(C == 0);
  ph1.add_con(A >= 0);
  ph1.add_con(3*A + 2*B > 0);
  ph1.add_con(3*A + 2*B > 0);

  Poly ph2(3, Topol::NNC);
  for (const auto& c : ph1.cons())
    ph2.add_con(c);

  Poly ph3(3, Topol::NNC, Spec_Elem::EMPTY);
  ph3.add_gen(point(A + B));
  for (const auto& g : ph2.gens())
    ph3.add_gen(g);

  return (ph1 == ph2 && ph2 == ph3);
}

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
END_MAIN
