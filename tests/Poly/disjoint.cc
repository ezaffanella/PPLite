/* PPLite: a lightweight library for convex polyhedra derived from PPL.
   Copyright (C) 2001-2010 Roberto <bagnara@cs.unipr.it>
   Copyright (C) 2010-2018 BUGSENG srl (http://bugseng.com)
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

Poly
half_strip(const Gen& p, const Linear_Expr& e, bool closed = true) {
  Linear_Expr e1(p.linear_expr());
  e1 += 3*Var(0);

  Gens gs;
  gs.push_back(p);
  gs.push_back(ray(e));
  if (p.is_point())
    gs.push_back(point(e1, p.divisor()));
  else {
    gs.push_back(closure_point(e1, p.divisor()));
    e1 -= Var(0);
    e1 += e.get(Var(1)) * p.divisor() * Var(1);
    gs.push_back(point(e1));
  }
  Poly ph(2, Spec_Elem::EMPTY, (closed ? Topol::CLOSED : Topol::NNC));
  ph.add_gens(gs);
  return ph;
}

bool
test01() {
  Var A(0);
  Var B(1);

  Poly ph1(2);
  ph1.add_con(B == 0);
  ph1.add_con(A >= 1);
  ph1.add_con(A <= 2);

  Poly ph2(2);
  ph2.add_con(A == 0);
  ph2.add_con(B >= 1);
  ph2.add_con(B <= 2);

  bool ok = ph1.is_disjoint_from(ph2);

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  return ok;
}

bool
test02() {
  Var A(0);
  Var B(1);

  Poly ph1 = half_strip(point(A + B), B);

  Poly ph2(2, Spec_Elem::EMPTY);
  ph2.add_gen(point(2*A + B));
  ph2.add_gen(point(4*A + 3*B));
  ph2.add_gen(ray(A - B));

  bool disjoint = ph1.is_disjoint_from(ph2);

  print_gens(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");

  return !disjoint;
}

bool
test03() {
  Var A(0);
  Var B(1);

  Poly ph1 = half_strip(point(A + B), B);
  Poly ph2 = half_strip(point(4*A + B), B);

  bool disjoint = ph1.is_disjoint_from(ph2);

  print_gens(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");

  return !disjoint;
}

bool
test04() {
  Var A(0);
  Var B(1);

  Poly ph1 = half_strip(point(A + B), B);
  Poly ph2 = half_strip(point(A + B), -B);

  bool disjoint = ph1.is_disjoint_from(ph2);

  print_gens(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");

  return !disjoint;
}

bool
test05() {
  Var A(0);
  Var B(1);

  Poly ph1 = half_strip(point(), B);

  Poly ph2(2, Spec_Elem::EMPTY);
  ph2.add_gen(point(2*A - 2*B));
  ph2.add_gen(point(-2*A + 2*B));
  ph2.add_gen(ray(-A - B));

  bool disjoint = ph1.is_disjoint_from(ph2);

  print_gens(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");

  return !disjoint;
}

bool
test06() {
  Var A(0);
  Var B(1);

  Poly ph1 = half_strip(point(A + B), B, false);

  Poly ph2(2, Spec_Elem::EMPTY, Topol::NNC);
  ph2.add_gen(point(3*A + B));
  ph2.add_gen(closure_point(2*A + B));
  ph2.add_gen(closure_point(4*A + 3*B));
  ph2.add_gen(ray(A - B));

  bool disjoint = ph1.is_disjoint_from(ph2);

  print_gens(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");

  return !disjoint;
}

bool
test07() {
  Var A(0);
  Var B(1);

  Poly ph1 = half_strip(point(A + B), B, false);
  Poly ph2 = half_strip(closure_point(4*A + B), B, false);

  bool disjoint = ph1.is_disjoint_from(ph2);

  print_gens(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");

  return disjoint;
}

bool
test08() {
  Var A(0);
  Var B(1);

  Poly ph1 = half_strip(point(A + B), B, false);
  Poly ph2 = half_strip(closure_point(A + B), -B, false);

  bool disjoint = ph1.is_disjoint_from(ph2);

  print_gens(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");

  return disjoint;
}

bool
test09() {
  Var A(0);
  Var B(1);

  Poly ph1 = half_strip(point(), B, false);

  Poly ph2(2, Spec_Elem::EMPTY, Topol::NNC);
  ph2.add_gen(point(-2*A - 2*B));
  ph2.add_gen(closure_point(2*A - 2*B));
  ph2.add_gen(closure_point(-2*A + 2*B));
  ph2.add_gen(ray(-A - B));

  bool disjoint = ph1.is_disjoint_from(ph2);

  print_gens(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");

  return disjoint;
}

bool
test10() {

  Var A(0); Var B(1);

  Poly x(2);
  x.add_con(A + B == 0);

  Poly y(2);
  y.add_con(A == 0);

  bool ok = (not x.is_disjoint_from(y));

  return ok;
}

bool
test11() {
  Poly ph1(2, Spec_Elem::EMPTY);
  Poly ph2(2, Spec_Elem::UNIVERSE);
  bool ok = ph1.is_disjoint_from(ph1)
    && ph1.is_disjoint_from(ph2)
    && ph2.is_disjoint_from(ph1)
    && not ph2.is_disjoint_from(ph2);
  return ok;
}

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
  DO_TEST(test03);
  DO_TEST(test04);
  DO_TEST(test05);
  DO_TEST(test06);
  DO_TEST(test07);
  DO_TEST(test08);
  DO_TEST(test09);
  DO_TEST(test10);
  DO_TEST(test11);
END_MAIN
