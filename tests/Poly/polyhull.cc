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
  gs1.push_back(ray(x));
  gs1.push_back(ray(y));

  Poly ph1(2, Spec_Elem::EMPTY);
  ph1.add_gens(gs1);

  print_gens(ph1, "*** ph1 ***");

  Gens gs2;
  gs2.push_back(point(-x + y));
  gs2.push_back(point(x + y));
  gs2.push_back(point(3*x));

  Poly ph2(2, Spec_Elem::EMPTY);
  ph2.add_gens(gs2);

  print_gens(ph2, "*** ph2 ***");

  Poly computed_result = ph1;

  computed_result.poly_hull_assign(ph2);

  Gens gs_known_result;
  gs_known_result.push_back(point());
  gs_known_result.push_back(point(-x + y));
  gs_known_result.push_back(ray(x));
  gs_known_result.push_back(ray(y));

  Poly known_result(2, Spec_Elem::EMPTY);
  known_result.add_gens(gs_known_result);

  print_gens(computed_result, "*** ph1.poly_hull_assign(ph2) ***");

  return computed_result == known_result;
}

bool
test02() {
  Var x(0);
  Var y(1);

  Poly ph1(2, Spec_Elem::EMPTY);

  Gens gs;
  gs.push_back(point());
  gs.push_back(ray(x + y));

  Poly ph2(2, Spec_Elem::EMPTY);
  ph2.add_gens(gs);

  print_gens(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");

  Poly computed_result1(ph1);

  computed_result1.poly_hull_assign(ph2);

  Poly known_result(ph2);

  bool ok = (computed_result1 == known_result);

  print_gens(computed_result1, "*** after poly_hull_assign ***");

  return ok;
}

bool
test03() {
  Var x(0);
  Var y(1);

  Poly ph1(2);
  ph1.add_con(x >= 0);
  ph1.add_con(y >= 0);
  ph1.add_con(x <= 2);
  ph1.add_con(y <= 2);

  Poly ph2(2);
  ph2.add_con(y >= 2);
  ph2.add_con(y <= 4);
  ph2.add_con(x >= 0);
  ph2.add_con(x <= 2);

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  ph1.poly_hull_assign(ph2);

  print_gens(ph1, "*** after poly_hull_assign ***");

  Poly known_result(2, Spec_Elem::EMPTY);
  known_result.add_gen(point());
  known_result.add_gen(point(2*x));
  known_result.add_gen(point(4*y));
  known_result.add_gen(point(2*x + 4*y));

  bool ok = (ph1 == known_result);

  return ok;
}

bool
aux_test04(Poly& ph1, const Poly& ph2,
           // Note intentional call-by-value!
           Poly known_result) {
  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  ph1.poly_hull_assign(ph2);

  print_gens(ph1, "*** after poly_hull_assign ***");

  return ph1 == known_result;
}

bool
test04() {
  Var x(0);
  Var y(1);

  Poly ph1_1(2);
  ph1_1.add_con(x >= 0);
  ph1_1.add_con(y >= 0);
  ph1_1.add_con(x <= 2);
  ph1_1.add_con(y <= 2);
  Poly ph1_2(ph1_1);

  Poly ph2_1(2);
  ph2_1.add_con(x+y <= 0);
  ph2_1.add_con(x+y >= 2);
  Poly ph2_2(ph2_1);
  Poly ph2_3(ph2_1);
  Poly ph2_4(ph2_1);

  bool ok = aux_test04(ph1_1, ph2_1, ph1_1)
    && aux_test04(ph2_2, ph1_2, ph1_2)
    && aux_test04(ph2_3, ph2_4, ph2_3);

  return ok;
}

bool
test05() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::EMPTY);
  ph1.add_gen(point(A));
  ph1.add_gen(ray(A));
  ph1.add_gen(ray(B));
  Poly ph2(2, Spec_Elem::EMPTY);

  print_gens(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");

  Poly known_result(ph1);
  ph1.poly_hull_assign(ph2);
  bool ok = (ph1 == known_result);
  print_gens(ph1, "*** after ph1.poly_hull_assign(ph2) ***");
  return ok;
}

bool
test06() {
  Poly ph1;
  Poly ph2;

  print_gens(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");

  Poly known_result(ph1);
  ph1.poly_hull_assign(ph2);
  bool ok = (ph1 == known_result);
  print_gens(ph1, "*** after ph1.poly_hull_assign(ph2) ***");
  return ok;
}

bool
test07() {
  Var A(0);
  Var B(1);

  Gens gs1;
  gs1.push_back(point());
  gs1.push_back(ray(A));
  gs1.push_back(point(B));
  Poly ph1(2, Spec_Elem::EMPTY);
  ph1.add_gens(gs1);
  ph1.minimize();

  Gens gs2;
  gs2.push_back(point());
  gs2.push_back(ray(B));
  Poly ph2(2, Spec_Elem::EMPTY);
  ph2.add_gens(gs2);

  print_cons(ph1, "*** ph1 ***");
  print_gens(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");

  ph1.poly_hull_assign(ph2);

  Poly known_result(2);
  known_result.add_con(A >= 0);
  known_result.add_con(B >= 0);
  bool ok = (ph1 == known_result);
  print_gens(ph1, "*** after ph1.poly_hull_assign(ph2) ***");
  return ok;
}

bool
test08() {
  Var A(0);
  Var B(1);

  Poly ph1(2);
  ph1.add_con(A - B >= 0);
  Poly ph2(2, Spec_Elem::EMPTY);

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  Poly known_result(ph1);
  ph1.poly_hull_assign(ph2);
  bool ok = (ph1 == known_result);
  print_cons(ph1, "*** after ph1.poly_hull_assign(ph2) ***");
  return ok;
}

bool
test09() {
  Var A(0);
  Var B(1);

  Poly ph1(2);
  ph1.minimize();
  ph1.add_con(A == B);
  Poly copy_ph1 = ph1;

  Poly ph2(2);
  ph2.minimize();
  ph2.add_con(A - B >= 1);
  Poly copy_ph2 = ph2;

  print_gens(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");

  ph1.poly_hull_assign(ph2);
  copy_ph1.poly_hull_assign(copy_ph2);
  bool ok = (ph1 == copy_ph1);
  print_gens(ph1, "*** after poly_hull_assign ***");
  print_gens(copy_ph1, "*** after poly_hull_assign ***");
  return ok;
}

bool
test10() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::EMPTY);
  ph1.add_gen(point());
  ph1.copy_cons();
  ph1.add_gen(line(A + B));
  Poly copy_ph1 = ph1;

  Poly ph2(2, Spec_Elem::EMPTY);
  ph2.add_gen(point());
  ph2.copy_cons();
  ph2.add_gen(ray(A));
  ph2.add_gen(ray(B));

  Poly copy_ph2 = ph2;

  print_gens(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");

  ph1.poly_hull_assign(ph2);
  copy_ph1.poly_hull_assign(copy_ph2);
  bool ok = (ph1 == copy_ph1);
  print_gens(ph1, "*** after poly_hull_assign ***");
  print_gens(copy_ph1, "*** after poly_hull_assign ***");
  return ok;
}

bool
test11() {
  Var B(1);
  Var C(2);

  Poly p(3);
  p.add_con(B >= 0);
  p.add_con(C >= 0);

  Poly q(3);
  q.add_con(C >= 0);

  print_cons(p, "*** p ***");
  print_cons(q, "*** q ***");

  p.poly_hull_assign(q);
  bool ok = (p == q);
  print_cons(p, "*** p.poly_hull_assign(q) ***");
  return ok;
}

bool
test12() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Topol::NNC);
  ph1.add_con(B == 0);
  ph1.add_con(A >= 0);
  ph1.add_con(B - A == 2);

  Poly ph2(2, Topol::NNC);
  ph2.add_con(B == 0);
  ph2.add_con(-A > 0);
  ph2.add_con(B - A == 0);

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  ph1.poly_hull_assign(ph2);

  bool ok = ph1.is_empty();

  print_cons(ph1, "*** ph1.poly_hull_assign(ph2) ***");

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
  DO_TEST(test12);
END_MAIN
