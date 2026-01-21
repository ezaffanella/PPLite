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
test00() {
  // Dummy test: this changes once and for all the widening implementation
  // to use the bhrz03 widening (rather than h79).
  set_widen_impl(Widen_Impl::BHRZ03);
  return true;
}

bool
test01() {
  Var A(0);
  Var B(1);

  Gens gs1;
  gs1.push_back(closure_point(A + B));
  gs1.push_back(closure_point(3*A + B));
  gs1.push_back(closure_point(A + 3*B));
  gs1.push_back(closure_point(3*A + 3*B));
  gs1.push_back(point(2*A + 2*B));
  Poly ph1(Topol::NNC, 2, Spec_Elem::EMPTY);
  ph1.add_gens(gs1);

  Gens gs2;
  gs2.push_back(closure_point());
  gs2.push_back(closure_point(4*A));
  gs2.push_back(closure_point(4*B));
  gs2.push_back(closure_point(4*A + 4*B));
  gs2.push_back(point(2*A));
  gs2.push_back(point(2*A + 4*B));
  gs2.push_back(point(2*B));
  gs2.push_back(point(4*A + 2*B));
  Poly ph2(Topol::NNC, 2, Spec_Elem::EMPTY);
  ph2.add_gens(gs2);

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  Poly known_result(ph2);

  // The "do not widen" technique applies.
  ph2.widening_assign(ph1);

  bool ok = (ph2 == known_result);

  print_cons(ph2, "*** ph2.widening_assign(ph1) ***");
  print_cons(known_result, "*** known_result ***");

  return ok;
}

bool
test02() {
  Var A(0);
  Var B(1);

  Gens gs1;
  gs1.push_back(point(A + B));
  gs1.push_back(ray(A));
  gs1.push_back(ray(B));
  Poly ph1(Topol::NNC, 2, Spec_Elem::EMPTY);
  ph1.add_gens(gs1);

  Gens gs2;
  gs2.push_back(closure_point());
  gs2.push_back(point(A));
  gs2.push_back(point(B));
  gs2.push_back(ray(A));
  gs2.push_back(ray(B));
  Poly ph2(Topol::NNC, 2, Spec_Elem::EMPTY);
  ph2.add_gens(gs2);

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  Poly known_result(Topol::NNC, 2);

  // Certificate convergence test fails: fall back to standard widening.
  ph2.widening_assign(ph1);

  bool ok = (ph2 == known_result);

  print_cons(ph2, "*** ph2.widening_assign(ph1) ***");
  print_cons(known_result, "*** known_result ***");

  return ok;
}

bool
test03() {
  Var A(0);
  Var B(1);

  Gens gs1;
  gs1.push_back(closure_point(A + B));
  gs1.push_back(point(2*A + 2*B));
  gs1.push_back(ray(A+B));
  gs1.push_back(point(2*A + B));
  gs1.push_back(ray(B));
  Poly ph1(Topol::NNC, 2, Spec_Elem::EMPTY);
  ph1.add_gens(gs1);

  Gens gs2;
  gs2.push_back(closure_point());
  gs2.push_back(point(A));
  gs2.push_back(point(B));
  gs2.push_back(ray(A));
  gs2.push_back(ray(B));
  Poly ph2(Topol::NNC, 2, Spec_Elem::EMPTY);
  ph2.add_gens(gs2);

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  Poly known_result(ph2);

  // Evolving rays.
  ph2.widening_assign(ph1);

  bool ok = (ph2 == known_result);

  print_cons(ph2, "*** ph2.widening_assign(ph1) ***");
  print_cons(known_result, "*** known_result ***");

  return ok;
}

bool
test04() {
  Var A(0);
  Var B(1);

  Gens gs1;
  gs1.push_back(point());
  gs1.push_back(point(A));
  gs1.push_back(point(B));
  gs1.push_back(point(A + B));
  Poly ph1(Topol::NNC, 2, Spec_Elem::EMPTY);
  ph1.add_gens(gs1);

  Gens gs2;
  gs2.push_back(point());
  gs2.push_back(point(A));
  gs2.push_back(point(-2*A + B));
  gs2.push_back(ray(A + B));
  gs2.push_back(ray(-A + B));
  Poly ph2(Topol::NNC, 2, Spec_Elem::EMPTY);
  ph2.add_gens(gs2);

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  Poly known_result(ph2);

  ph2.widening_assign(ph1);

  bool ok = (ph2 == known_result);

  print_cons(ph2, "*** ph2.widening_assign(ph1) ***");
  print_cons(known_result, "*** known_result ***");

  return ok;
}

bool
test05() {
  Var A(0);
  Var B(1);

  Cons cs1;
  cs1.push_back(B >= 1);
  cs1.push_back(A >= 0);
  cs1.push_back(A < 2);
  cs1.push_back(A + B > 1);
  Poly ph1(Topol::NNC, 2);
  ph1.add_cons(cs1);

  Cons cs2;
  cs2.push_back(B > 0);
  cs2.push_back(A >= 0);
  cs2.push_back(A <= 2);
  Poly ph2(Topol::NNC, 2);
  ph2.add_cons(cs2);

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  Poly known_result(ph2);

  ph2.widening_assign(ph1);

  bool ok = (ph2 == known_result);

  print_cons(ph2, "*** ph2.widening_assign(ph1) ***");
  print_cons(known_result, "*** known_result ***");

  return ok;
}

bool
test06() {
  Var A(0);
  Var B(1);

  Cons cs1;
  cs1.push_back(B >= 1);
  cs1.push_back(A >= 0);
  cs1.push_back(A < 2);
  cs1.push_back(A + B > 1);
  cs1.push_back(B <= 3);
  Poly ph1(Topol::NNC, 2);
  ph1.add_cons(cs1);

  Cons cs2;
  cs2.push_back(B > 0);
  cs2.push_back(A >= 0);
  cs2.push_back(A <= 2);
  cs2.push_back(B <= 3);
  Poly ph2(Topol::NNC, 2);
  ph2.add_cons(cs2);

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  Poly known_result(ph2);

  ph2.widening_assign(ph1);

  bool ok = (ph2 == known_result);

  print_cons(ph2, "*** ph2.widening_assign(ph1) ***");
  print_cons(known_result, "*** known_result ***");

  return ok;
}

bool
test08() {
  Var A(0);
  Var B(1);

  Poly ph1(Topol::NNC, 2, Spec_Elem::EMPTY);
  ph1.add_gen(point());
  ph1.add_gen(closure_point(A));
  ph1.add_gen(closure_point(B));
  ph1.add_gen(closure_point(A + B));
  ph1.minimize();

  Poly ph2(Topol::NNC, 2, Spec_Elem::EMPTY);
  ph2.add_gen(point());
  ph2.add_gen(point(A));
  ph2.add_gen(point(B));
  ph2.add_gen(closure_point(2*A + 2*B));
  ph2.minimize();

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  Poly known_result(ph2);

  ph2.widening_assign(ph1);

  bool ok = (ph2 == known_result);

  print_cons(ph2, "*** ph2.widening_assign(ph1) ***");
  print_cons(known_result, "*** known_result ***");

  return ok;
}

bool
test09() {
  Var A(0);
  Var B(1);

  Poly ph1(Topol::NNC, 2, Spec_Elem::EMPTY);
  ph1.add_gen(point());
  ph1.add_gen(point(A));
  ph1.add_gen(point(B));
  ph1.add_gen(closure_point(A + B));
  ph1.minimize();

  Poly ph2(Topol::NNC, 2, Spec_Elem::EMPTY);
  ph2.add_gen(point());
  ph2.add_gen(point(A));
  ph2.add_gen(point(B));
  ph2.add_gen(closure_point(2*A + 2*B));
  ph2.minimize();

  Poly ph3(Topol::NNC, 2, Spec_Elem::EMPTY);
  ph3.add_gen(point());
  ph3.add_gen(point(A));
  ph3.add_gen(point(B));
  ph3.add_gen(ray(2*A + B));
  ph3.add_gen(ray(A + 2*B));
  ph3.minimize();

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  Poly known_result(ph3);

  ph2.widening_assign(ph1);

  bool ok = (ph2 == known_result);

  print_cons(ph2, "*** ph2.widening_assign(ph1) ***");
  print_cons(known_result, "*** known_result ***");

  return ok;
}

BEGIN_MAIN
  DO_TEST(test00);
  DO_TEST(test01);
  DO_TEST(test02);
  DO_TEST(test03);
  DO_TEST(test04);
  DO_TEST(test05);
  DO_TEST(test06);
  DO_TEST(test08);
  DO_TEST(test09);
END_MAIN
