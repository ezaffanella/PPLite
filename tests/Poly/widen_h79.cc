/* PPLite: a lightweight library for convex polyhedra derived from PPL.
   Copyright (C) 2001-2010 Roberto <bagnara@cs.unipr.it>
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

#include "pplite_test.hh"

bool
test01() {
  Var x(0);
  Var y(1);

  Poly ph1(2);
  ph1.add_con(x >= 0);
  ph1.add_con(y >= 0);
  ph1.add_con(x <= 2);
  ph1.add_con(y <= 2);

  Poly ph2(2);
  // Add inconsistent constraints to ph2.
  ph2.add_con(x+y <= 0);
  ph2.add_con(x+y >= 2);

  Poly ph1_1(ph1);
  print_cons(ph1_1, "*** ph1_1 ***");
  Poly ph2_1(ph2);
  print_cons(ph2_1, "*** ph2_1 ***");
  ph1_1.widening_assign(ph2_1);
  print_gens(ph1_1, "*** after widening_assign ***");
  Poly ph1_2(ph1);
  bool ok = (ph1_1 == ph1_2);

  Poly ph2_2(ph2);
  print_cons(ph2_2, "*** ph2_2 ***");
  Poly ph2_3(ph2);
  print_cons(ph2_3, "*** ph2_3 ***");
  ph2_2.widening_assign(ph2_3);
  print_gens(ph2_2, "*** after widening_assign ***");
  Poly ph2_4(ph2);
  ok = ok && (ph2_2 == ph2_4);

  return ok;
}

bool
test02() {
  Var A(0);
  Var B(1);

  Poly ph1(2);
  ph1.add_con(A >= 2);
  ph1.add_con(B >= 0);

  Poly ph2(2);
  ph2.add_con(A >= 0);
  ph2.add_con(B >= 0);
  ph2.add_con(A-B >= 2);

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  ph1.widening_assign(ph2);

  Poly known_result(2);
  known_result.add_con(B >= 0);

  bool ok = (ph1 == known_result);

  print_cons(ph1, "*** after widening_assign ***");

  return ok;
}

bool
test03() {
  Var A(0);
  Var B(1);

  Gens gs1;
  gs1.push_back(point());
  gs1.push_back(ray(A));
  gs1.push_back(ray(A + B));
  Poly ph1(2, Spec_Elem::EMPTY);
  ph1.add_gens(gs1);

  Gens gs2;
  gs2.push_back(point());
  gs2.push_back(ray(A));
  gs2.push_back(ray(A + 2*B));
  Poly ph2(2, Spec_Elem::EMPTY);
  ph2.add_gens(gs2);

  print_gens(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");

  ph2.widening_assign(ph1);

  Poly known_result(2);
  known_result.add_con(B >= 0);

  bool ok = (ph2 == known_result);

  print_cons(ph2, "*** after ph2.widening_assign(ph1) ***");

  return ok;
}

bool
test04() {
  Poly ph1;
  Poly ph2(0, Spec_Elem::EMPTY);

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  ph1.widening_assign(ph2);

  Poly known_result;
  known_result = ph1;

  bool ok = (ph1 == known_result);

  print_cons(ph1, "*** after ph1.widening_assign(ph2) ***");

  return ok;
}

bool
test05() {
  Var A(0);
  Var B(1);

  Poly ph1(2);
  ph1.add_con(-A - 2*B >= -6);
  ph1.add_con(B >= 0);
  ph1.add_con(A - 2*B >= 2);

  Poly ph2(2);
  ph2.add_con(-A - 2*B >= -10);
  ph2.add_con(B >= 0);
  ph2.add_con(A - 2*B >= 2);

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  ph2.widening_assign(ph1);

  Poly known_result(2);
  known_result.add_con(B >= 0);
  known_result.add_con(A - 2*B >= 2);

  bool ok = (ph2 == known_result);

  print_cons(ph2, "*** after ph2.widening_assign(ph1) ***");

  return ok;
}

bool
test06() {
  Var A(0);
  Var B(1);

  Poly ph1(2);
  ph1.add_con(-A + B >= 2);
  ph1.add_con(A >= 0);

  Poly ph2(2);
  ph2.add_con(-A + B >= 3);
  ph2.add_con(A - B >= -8);
  ph2.add_con(A >= 1);

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  ph1.widening_assign(ph2);

  bool ok = ph1.is_universe();

  print_cons(ph1, "*** after widening_assign ***");

  return ok;
}

bool
test07() {
  Var A(0);
  Var B(1);

  Poly ph1(2);
  ph1.add_con(A >= 0);
  ph1.add_con(B >= 0);
  ph1.add_gen(line(A));

  Poly ph2(2);
  ph2.add_con(A >= 0);
  ph2.add_con(B >= 0);

  Poly ph1_copy = ph1;
  print_cons(ph1_copy, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  ph1.widening_assign(ph2);

  Poly known_result(2);
  known_result.add_con(B >= 0);

  bool ok = (ph1 == known_result);

  print_cons(ph1, "*** after widening_assign ***");

  return ok;
}

bool
test08() {
  Var A(0);

  Poly ph1(1, Topol::NNC);
  ph1.add_con(A > 0);
  ph1.add_con(A < 2);

  Poly ph4(1, Topol::NNC);
  ph4.add_con(4*A >= 1);
  ph4.add_con(4*A <= 3);

  Poly ph = ph4;
  ph.intersection_assign(ph1);

  print_cons(ph4, "*** ph4 ***");
  print_cons(ph, "*** ph ***");

  Poly known_result(ph4);

  ph.widening_assign(ph4);

  bool ok = (ph == known_result);

  print_cons(ph, "*** after widening_assign ***");

  return ok;
}

bool
test09() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Topol::NNC);
  ph1.add_con(A >= 0);
  ph1.add_con(A <= 2);
  ph1.add_con(B >= 0);
  ph1.add_con(B <= 2);
  ph1.add_con(A + B > 0);

  Poly ph2(2, Topol::NNC);
  ph2.add_con(A >= 0);
  ph2.add_con(B >= 0);
  ph2.add_con(3*A + 2*B > 0);

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  Poly known_result(ph2);

  ph2.widening_assign(ph1);

  bool ok = (ph2 == known_result);

  print_cons(ph2, "*** ph2.widening_assign(ph1) ***");

  return ok;
}

bool
test10() {
  Var A(0);
  Var B(1);
  Poly ph1(2);
  ph1.add_con(A == 0);
  ph1.add_con(B >= 0);
  Poly ph2(2);
  ph2.add_con(A == 0);
  ph2.add_con(B == 0);

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  Poly known_result(ph1);
  ph1.widening_assign(ph2);
  bool ok = (ph1 == known_result);
  print_cons(ph1, "*** ph1.widening_assign(ph2) ***");
  return ok;
}

bool
test11() {
  Var A(0);
  Var B(1);
  Var C(2);
  Var D(3);

  Poly ph1(4, Topol::NNC);
  ph1.add_con(A <= 19);
  ph1.add_con(B >= 0);
  ph1.add_con(B <= 1);
  ph1.add_con(D < 0);
  ph1.minimize();

  Poly ph2(4, Topol::NNC);
  ph2.add_con(A == 1);
  ph2.add_con(B == 0);
  ph2.add_con(10*C - D == 0);
  ph2.add_con(C < 0);
  ph2.minimize();

  print_cons(ph1, "=== ph1 ===");
  print_cons(ph2, "=== ph2 ===");

  if (not (ph1.check_inv() && ph2.check_inv() && ph1.contains(ph2)))
    return false;

  Poly ph(ph1);
  ph.widening_assign(ph2);
  ph.minimize();

  print_cons(ph, "=== ph1.widening_assign(ph2) ===");
  if (not (ph.check_inv() && ph.contains(ph1) && ph.contains(ph2)))
    return false;

  // This would be the proper result of h79; we are more precise
  // due to the check for affine dimension.
  // Poly kres(4, Topol::NNC);
  // kres.add_con(B >= 0);
  // kres.add_con(D < 0);
  // kres.minimize();

  Poly kres = ph1;
  return (ph == kres);
}

bool
test12() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Topol::NNC);
  ph1.add_con(B < 0);
  ph1.minimize();

  Poly ph2(2, Topol::NNC);
  ph2.add_con(2*A - B == 0);
  ph2.add_con(A < 0);
  ph2.minimize();

  print_cons(ph1, "=== ph1 ===");
  print_cons(ph2, "=== ph2 ===");

  if (not (ph1.check_inv() && ph2.check_inv() && ph1.contains(ph2)))
    return false;

  Poly ph(ph1);
  ph.widening_assign(ph2);
  ph.minimize();

  print_cons(ph, "=== ph1.widening_assign(ph2) ===");
  if (not (ph.check_inv() && ph.contains(ph1) && ph.contains(ph2)))
    return false;

  Poly kres(2, Topol::NNC);
  kres.add_con(B < 0);
  kres.minimize();

  return (ph == kres);
}

bool
test13() {
  auto old_spec = get_widen_spec();
  set_widen_spec(Widen_Spec::SAFE);

  Var A(0);

  Poly ph1(2, Topol::NNC);
  ph1.add_con(A == 2);
  ph1.minimize();

  Poly ph2(2, Topol::NNC);
  ph2.add_con(A >= 1);
  ph2.add_con(A <= 2);
  ph2.minimize();

  print_cons(ph1, "=== ph1 ===");
  print_cons(ph2, "=== ph2 ===");

  // Safe widening does not have the inclusion precondition.
  set_widen_impl(Widen_Impl::H79);
  Poly ph_safe = ph1;
  Poly kres = ph2;
  ph_safe.widening_assign(ph2);

  print_cons(ph_safe, "=== ph1.widening_assign(ph2) safe ===");

  set_widen_spec(old_spec);

  return (ph_safe == kres);
}

BEGIN_MAIN
  set_widen_spec(Widen_Spec::RISKY);
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
  DO_TEST(test13);
END_MAIN
