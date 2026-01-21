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
  Var x(0);
  Var y(1);

  Gens gs1;
  gs1.push_back(point());
  gs1.push_back(point(4*x));
  gs1.push_back(point(2*x + 2*y));

  Poly ph1(2, Spec_Elem::EMPTY);
  ph1.add_gens(gs1);

  print_gens(ph1, "*** ph1 ***");

  Gens gs2;
  gs2.push_back(point(0*x + 3*y));
  gs2.push_back(point(4*x + 3*y));
  gs2.push_back(point(2*x + 1*y));

  Poly ph2(2, Spec_Elem::EMPTY);
  ph2.add_gens(gs2);

  print_gens(ph2, "*** ph2 ***");

  Poly computed_result = ph1;

  computed_result.poly_difference_assign(ph2);

  Gens gs_known_result;
  gs_known_result.push_back(point());
  gs_known_result.push_back(point(3*x + 3*y, 2));
  gs_known_result.push_back(point(4*x));
  gs_known_result.push_back(point(5*x + 3*y, 2));

  Poly known_result(2, Spec_Elem::EMPTY);
  known_result.add_gens(gs_known_result);

  Poly ph3(2);
  ph3.add_con(2*y >= 3);

  known_result.poly_difference_assign(ph3);

  bool ok = (computed_result == known_result);

  print_gens(computed_result, "*** after poly_difference_assign ***");
  print_gens(known_result, "*** known_result ***");

  return ok;
}

bool
test02() {
  Poly ph1(0);
  Poly ph2(0);

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  Poly known_result(0);
  known_result.add_con(Linear_Expr() <= -4);

  ph1.poly_difference_assign(ph2);

  bool ok = (ph1 == known_result);

  print_cons(ph1, "*** after ph1.poly_difference_assign(ph2) ***");

  return ok;
}

bool
test03() {
  Var A(0);
  Var B(1);

  Poly ph1(2);
  ph1.add_con(A >= 0);
  ph1.add_con(A <= 2);
  ph1.add_con(B == 0);

  Poly ph2(2);
  ph2.add_con(A >= 0);
  ph2.add_con(A <= 2);
  ph2.add_con(B >= 0);
  ph2.add_con(B <= 2);

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  Poly known_result(2, Spec_Elem::EMPTY);

  ph1.poly_difference_assign(ph2);

  bool ok = (ph1 == known_result);

  print_cons(ph1, "*** after ph1.poly_difference_assign(ph2) ***");

  return ok;
}

bool
test04() {
  Var A(0);

  Poly ph1(1);
  ph1.add_con(A >= 0);
  ph1.add_con(A <= 7);
  Poly ph2(1);
  ph2.add_con(A == 5);

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  Poly known_result(ph1);

  ph1.poly_difference_assign(ph2);

  bool ok = (ph1 == known_result);

  print_cons(ph1, "*** after ph1.poly_difference_assign(ph2) ***");

  return ok;
}

bool
test05() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::EMPTY);
  Poly ph2(2);
  ph2.add_con(A == B);

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  Poly known_result(ph1);

  ph1.poly_difference_assign(ph2);

  bool ok = (ph1 == known_result);

  print_cons(ph1, "*** after ph1.poly_difference_assign(ph2) ***");

  return ok;
}

bool
test06() {
  Var A(0);
  Var B(1);

  Poly ph1(2);
  ph1.add_con(A + B == 0);
  Poly ph2(2, Spec_Elem::EMPTY);

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  Poly known_result(ph1);

  ph1.poly_difference_assign(ph2);

  bool ok = (ph1 == known_result);

  print_cons(ph1, "*** after ph1.poly_difference_assign(ph2) ***");

  return ok;
}

bool
test07() {
  Var x(0);
  Var y(1);

  Gens gs1;
  gs1.push_back(point());
  gs1.push_back(point(3*y));
  gs1.push_back(point(3*x));
  gs1.push_back(point(3*x + 3*y));

  Poly ph1(2, Spec_Elem::EMPTY, Topol::NNC);
  ph1.add_gens(gs1);

  print_cons(ph1, "*** ph1 ***");

  Cons cs2;
  cs2.push_back(x == 0);

  Poly ph2(2, Topol::NNC);
  ph2.add_cons(cs2);

  print_cons(ph2, "*** ph2 ***");

  Poly computed_result = ph1;
  computed_result.poly_difference_assign(ph2);

  Gens gs_known_result;
  gs_known_result.push_back(closure_point());
  gs_known_result.push_back(closure_point(3*y));
  gs_known_result.push_back(point(3*x));
  gs_known_result.push_back(point(3*x + 3*y));

  Poly known_result(2, Spec_Elem::EMPTY, Topol::NNC);
  known_result.add_gens(gs_known_result);

  bool ok = (computed_result == known_result);

  print_gens(computed_result, "*** after difference_assign ***");
  print_gens(known_result, "*** known_result ***");

  return ok;
}

bool
test08() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Topol::NNC);
  ph1.add_con(A >= 0);
  ph1.add_con(B >= 0);

  Poly ph2(2, Topol::NNC);
  ph2.add_con(A > 2);
  ph2.add_con(B >= 0);

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  ph1.poly_difference_assign(ph2);

  Poly known_result(2, Topol::NNC);
  known_result.add_con(A >= 0);
  known_result.add_con(A <= 2);
  known_result.add_con(B >= 0);

  bool ok = (ph1 == known_result);

  print_cons(ph1, "*** after ph1.poly_difference_assign(ph2) ***");

  return ok;
}

bool
test09() {
  Var x(0);
  Var y(1);

  Cons cs1;
  cs1.push_back(x >= 0);
  cs1.push_back(y >= 0);
  cs1.push_back(x <= 4);
  cs1.push_back(x - 2*y <= 2);

  Poly ph1(2);
  ph1.add_cons(cs1);

  print_cons(ph1, "*** ph1 ***");

  Cons cs2;
  cs2.push_back(x >= 0);
  cs2.push_back(y >= 0);
  cs2.push_back(x <= 4);
  cs2.push_back(y <= 5);
  cs2.push_back(x - 2*y <= 2);
  cs2.push_back(x + y <= 7);

  Poly ph2(2);
  ph2.add_cons(cs2);

  print_cons(ph2, "*** ph2 ***");

  Poly computed_result = ph1;
  computed_result.poly_difference_assign(ph2);

  Gens gs_known_result;
  gs_known_result.push_back(point(0*x + 5*y));
  gs_known_result.push_back(point(4*x + 3*y));
  gs_known_result.push_back(ray(0*x + 1*y));

  Poly known_result(2, Spec_Elem::EMPTY);
  known_result.add_gens(gs_known_result);

  bool ok = (computed_result == known_result);

  print_cons(computed_result, "*** after difference_assign ***");
  print_cons(known_result, "*** known_result ***");

  return ok;
}

bool
test10() {
  Var A(0);

  Poly ph1(1);

  Poly ph2(1);
  ph2.add_con(A >= 1);
  ph2.add_con(A <= 0);

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  ph1.poly_difference_assign(ph2);

  Poly known_result(1);

  bool ok = (ph1 == known_result);

  print_cons(ph1, "*** after ph1.poly_difference_assign(ph2) ***");

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
END_MAIN
