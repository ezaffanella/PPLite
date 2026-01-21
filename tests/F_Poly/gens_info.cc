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
test01() {
  Var A(0), B(1), C(2);
  Cons cs;
  cs.push_back(A > 0);
  cs.push_back(A < 2);
  cs.push_back(B > 0);
  cs.push_back(B < 2);
  cs.push_back(C > 0);
  cs.push_back(C < 2);
  Gens gs;
  gs.push_back(point(A));
  gs.push_back(point(B));
  gs.push_back(point(A + 2*C));
  gs.push_back(point(B + 2*C));

  Poly ph1(3, Topol::NNC);
  ph1.add_cons(cs);
  ph1.add_gens(gs);
  auto gi1 = ph1.gens_info();

  F_Poly ph2(3, Topol::NNC);
  ph2.add_cons(cs);
  ph2.add_gens(gs);
  auto gi2 = ph2.gens_info();

  return gi1 == gi2;
}

bool
test02() {
  Var A(0), B(1), C(2);
  Cons cs;
  cs.push_back(A > 0);
  cs.push_back(A < 2);
  cs.push_back(B > 0);
  cs.push_back(B < 2);
  Gens gs;
  gs.push_back(point(A));
  gs.push_back(point(B));

  Poly ph1(2, Topol::NNC);
  ph1.add_cons(cs);
  ph1.add_gens(gs);
  Poly tmp1(1, Topol::NNC);
  tmp1.add_con(A >= 0);
  tmp1.add_con(A <= 1);
  ph1.concatenate_assign(tmp1);
  auto gi1 = ph1.gens_info();

  F_Poly ph2(2, Topol::NNC);
  ph2.add_cons(cs);
  ph2.add_gens(gs);
  F_Poly tmp2(1, Topol::NNC);
  tmp2.add_con(A >= 0);
  tmp2.add_con(A <= 1);
  ph2.concatenate_assign(tmp2);
  auto gi2 = ph2.gens_info();

  return gi1 == gi2;
}

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
END_MAIN
