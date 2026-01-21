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
  Var C(2);

  Gens gs1;
  gs1.push_back(point());
  gs1.push_back(point(A));
  gs1.push_back(point(B));
  gs1.push_back(point(A + B));
  gs1.push_back(point(C));
  gs1.push_back(point(A + C));
  gs1.push_back(point(B + C));
  gs1.push_back(point(A + B + C));
  Poly ph1(3, Spec_Elem::EMPTY);
  ph1.add_gens(gs1);

  Poly ph1_copy(ph1);

  Gens gs2;
  gs2.push_back(point());
  gs2.push_back(point(A));
  gs2.push_back(point(B));
  gs2.push_back(point(A + B));
  Poly ph2(3, Spec_Elem::EMPTY);
  ph2.add_gens(gs2);

  print_cons(ph1, "*** ph1 ***");
  print_gens(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");
  print_gens(ph2, "*** ph2 ***");

  ph1.widening_assign(ph2, Widen_Impl::BHRZ03);

  bool ok = (ph1 == ph1_copy);

  print_cons(ph1, "*** after widening_assign ***");
  print_cons(ph1_copy, "*** expected result ***");

  return ok;
}

bool
test02() {
  Var A(0);
  Var B(1);

  Poly ph1(2);
  ph1.add_con(A >= 0);

  Poly ph1_copy(ph1);

  Poly ph2(2);
  ph2.add_con(A >= 0);
  ph2.add_con(B >= 0);

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  ph1.widening_assign(ph2, Widen_Impl::BHRZ03);

  bool ok = (ph1 == ph1_copy);

  print_cons(ph1, "*** after widening_assign ***");

  return ok;
}

bool
test03() {
  Var A(0);
  Var B(1);

  Gens gs1;
  gs1.push_back(point());
  gs1.push_back(point(2*A));
  gs1.push_back(point(2*B));

  Poly ph1(2, Spec_Elem::EMPTY);
  ph1.add_gens(gs1);
  Poly ph1_copy(ph1);

  Gens gs2;
  gs2.push_back(point());
  gs2.push_back(point(A));
  gs2.push_back(point(B));
  gs2.push_back(point(A + B));
  Poly ph2(2, Spec_Elem::EMPTY);
  ph2.add_gens(gs2);

  print_gens(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");

  ph1.widening_assign(ph2, Widen_Impl::BHRZ03);

  bool ok = (ph1 == ph1_copy);

  print_cons(ph1, "*** after widening_assign ***");

  return ok;
}

bool
test04() {
  Var A(0);
  Var B(1);

  Gens gs1;
  gs1.push_back(point(2*A));
  gs1.push_back(closure_point(A+B));
  gs1.push_back(closure_point(3*A+B));
  Poly ph1(Topol::NNC, 2, Spec_Elem::EMPTY);
  ph1.add_gens(gs1);

  Gens gs2;
  gs2.push_back(point(2*A));
  gs2.push_back(closure_point(B));
  gs2.push_back(closure_point(4*A+B));
  Poly ph2(Topol::NNC, 2, Spec_Elem::EMPTY);
  ph2.add_gens(gs2);

  print_cons(ph1, "*** ph1 ***");
  print_gens(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");
  print_gens(ph2, "*** ph2 ***");

  ph2.widening_assign(ph1, Widen_Impl::BHRZ03);

  Poly known_result(Topol::NNC, 2);
  known_result.add_con(B >= 0);
  known_result.add_con(B < 1);

  bool ok = (ph2 == known_result);

  print_cons(ph2, "*** after widening_assign ***");

  return ok;
}

bool
test05() {
  Var A(0);
  Var B(1);

  Poly ph1(2);
  ph1.add_con(A - B >= 0);
  ph1.add_con(A + B <= 2);
  ph1.add_con(B >= 0);

  Poly ph2(2);
  ph2.add_con(2*A - B >= 0);
  ph2.add_con(B >= 0);
  ph2.add_con(2*A + B <= 4);

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  ph2.widening_assign(ph1, Widen_Impl::BHRZ03);

  Poly known_result(2);
  known_result.add_con(2*A + B >= 0);
  known_result.add_con(2*A - B <= 4);
  known_result.add_con(B >= 0);

  bool ok = (ph2 == known_result);

  print_cons(ph2, "*** after widening_assign ***");

  return ok;
}

bool
test06() {
  Var A(0);
  Var B(1);

  Poly ph1(2);
  ph1.add_con(A >= 0);
  ph1.add_con(A - B <= 0);

  Poly ph2(2);
  ph2.add_con(A >= 0);
  ph2.add_con(2*A - B <= 0);

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  ph1.widening_assign(ph2, Widen_Impl::BHRZ03);

  Poly known_result(2);
  known_result.add_con(B >= 0);
  known_result.add_con(A >= 0);

  bool ok = (ph1 == known_result);

  print_cons(ph1, "*** after widening_assign ***");

  return ok;
}

bool
test07() {
  Var A(0);
  Var B(1);
  Var C(2);

  Gens gs1;
  gs1.push_back(point());
  gs1.push_back(ray(A));
  gs1.push_back(ray(B));
  gs1.push_back(ray(A + 4*B + 2*C));
  Poly ph1(3, Spec_Elem::EMPTY);
  ph1.add_gens(gs1);

  Gens gs2;
  gs2.push_back(point());
  gs2.push_back(ray(A));
  gs2.push_back(ray(B));
  gs2.push_back(ray(A + 2*B + 4*C));
  Poly ph2(3, Spec_Elem::EMPTY);
  ph2.add_gens(gs2);

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  ph2.widening_assign(ph1, Widen_Impl::BHRZ03);

  Poly known_result(3);
  known_result.add_con(A >= 0);
  known_result.add_con(B >= 0);
  known_result.add_con(C >= 0);

  bool ok = (ph2 == known_result);

  print_cons(ph2, "*** after widening_assign ***");

  return ok;
}

bool
test08() {
  Var A(0);
  Var B(1);
  Var C(2);

  Gens gs1;
  gs1.push_back(point());
  gs1.push_back(line(A));
  gs1.push_back(ray(B));
  gs1.push_back(ray(A + B + C));
  Poly ph1(3, Spec_Elem::EMPTY);
  ph1.add_gens(gs1);

  Gens gs2;
  gs2.push_back(point());
  gs2.push_back(line(A));
  gs2.push_back(ray(B));
  gs2.push_back(ray(A + B + 2*C));
  Poly ph2(3, Spec_Elem::EMPTY);
  ph2.add_gens(gs2);

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  ph2.widening_assign(ph1, Widen_Impl::BHRZ03);

  Poly known_result(3);
  known_result.add_con(B >= 0);
  known_result.add_con(C >= 0);

  bool ok = (ph2 == known_result);

  print_cons(ph2, "*** after widening_assign ***");

  return ok;
}

Gens
aux1_test09() {
  Var A(0);
  Var B(1);
  Var C(2);
  Var D(3);
  Var E(4);

  Gens gs;
  gs.push_back(point());
  gs.push_back(ray(C));
  gs.push_back(ray(D));
  gs.push_back(ray(E));
  gs.push_back(ray(A + D));
  gs.push_back(ray(A + B + E));
  return gs;
}

Poly
aux2_test09(unsigned n) {
  Var A(0);
  Var B(1);
  Var C(2);
  Var D(3);
  Var E(4);

  Poly ph(5, Spec_Elem::EMPTY);
  ph.add_gens(aux1_test09());
  n += 2;
  ph.add_gen(ray(A + (n-1)*B + E));
  if (n % 2 == 0) {
    // Even.
    unsigned m = n / 2;
    ph.add_gen(ray(m*B + E));
    ph.add_gen(ray(A + (m-1)*B + D));
  }
  else {
    // Odd.
    ph.add_gen(ray(n*B + 2*E));
    ph.add_gen(ray(2*A + (n-2)*B + 2*D));
  }
  return ph;
}

bool
test09() {
  // Chain condition for widenings:
  // for each increasing chain of descriptions p_0, p_1, ..., p_i, ...,
  // the sequence q_0, q_1, ..., q_i, ... defined by q_0 = p_0 and
  // for each i >= 1, q_i = q_{i-1} \nabla p_i is ultimately stationary.
  Var A(0);
  Var B(1);
  Var C(2);
  Var D(3);
  Var E(4);

  // Initialization: set q_0.
  Poly q_i_minus_1 = aux2_test09(0);

  for (unsigned i = 1; i <= 100; ++i) {
    print_gens(q_i_minus_1, "*** Result of the previous iteration ***");

    Poly p_i = aux2_test09(i);
    print_gens(p_i, "*** New stuff ***");

    Poly q_i = q_i_minus_1;
    q_i.poly_hull_assign(p_i);
    print_gens(q_i, "*** Poly-hull of previous with new ***");

    q_i.widening_assign(q_i_minus_1, Widen_Impl::BHRZ03);
    print_gens(q_i, "*** Result of widening poly-hull with new ***");

    if (q_i == q_i_minus_1) {

      Poly known_result(5);
      known_result.add_con(A >= 0);
      known_result.add_con(B >= 0);
      known_result.add_con(C >= 0);
      known_result.add_con(D >= 0);
      known_result.add_con(E >= 0);
      known_result.add_con(-A + B + D >= 0);

      bool ok = (q_i == known_result);

      print_cons(q_i, "*** The cons of the fix point ***");
      print_gens(q_i, "*** The gens of the fix point ***");

      return ok;
    }
    q_i_minus_1 = q_i;
  }
  return false;
}

bool
test10() {
  Poly ph1(0);
  Poly ph2(0, Spec_Elem::EMPTY);

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  ph1.widening_assign(ph2, Widen_Impl::BHRZ03);

  Poly known_result(0);

  bool ok = (ph1 == known_result);

  print_cons(ph1, "*** after ph1.widening_assign(ph2) ***");

  return ok;
}

bool
test11() {
  Var A(0);
  Var B(1);

  Poly ph1(2);
  ph1.add_con(B >= 0);

  Poly ph2(2);
  ph2.add_con(A >= 2);
  ph2.add_con(A <= 0);

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  Poly known_result(ph1);

  ph1.widening_assign(ph2, Widen_Impl::BHRZ03);

  bool ok = (ph1 == known_result);

  print_cons(ph1, "*** after ph1.widening_assign(ph2) ***");

  return ok;
}

bool
test12() {
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

  ph2.widening_assign(ph1, Widen_Impl::BHRZ03);

  Poly known_result(2);
  known_result.add_con(A >= 0);
  known_result.add_con(B >= 0);

  bool ok = (ph2 == known_result);

  print_gens(ph2, "*** after ph2.widening_assign(ph1) ***");

  return ok;
}

bool
test13() {
  Var A(0);
  Var B(1);

  Gens gs1;
  gs1.push_back(point());
  gs1.push_back(ray(-A));
  gs1.push_back(ray(-A + B));
  Poly ph1(2, Spec_Elem::EMPTY);
  ph1.add_gens(gs1);

  Gens gs2;
  gs2.push_back(point());
  gs2.push_back(ray(-A));
  gs2.push_back(ray(-A + 2*B));
  Poly ph2(2, Spec_Elem::EMPTY);
  ph2.add_gens(gs2);

  print_gens(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");

  ph2.widening_assign(ph1, Widen_Impl::BHRZ03);

  Poly known_result(2);
  known_result.add_con(A <= 0);
  known_result.add_con(B >= 0);

  bool ok = (ph2 == known_result);

  print_gens(ph2, "*** after ph2.widening_assign(ph1) ***");

  return ok;
}

bool
test14() {
  Var A(0);
  Var B(1);

  Gens gs1;
  gs1.push_back(point());
  gs1.push_back(ray(-A));
  gs1.push_back(ray(-A - B));
  Poly ph1(2, Spec_Elem::EMPTY);
  ph1.add_gens(gs1);

  Gens gs2;
  gs2.push_back(point());
  gs2.push_back(ray(-A));
  gs2.push_back(ray(-A - 2*B));
  Poly ph2(2, Spec_Elem::EMPTY);
  ph2.add_gens(gs2);

  print_gens(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");

  ph2.widening_assign(ph1, Widen_Impl::BHRZ03);

  Poly known_result(2);
  known_result.add_con(A <= 0);
  known_result.add_con(B <= 0);

  bool ok = (ph2 == known_result);

  print_gens(ph2, "*** after ph2.widening_assign(ph1) ***");

  return ok;
}

bool
test15() {
  Var A(0);
  Var B(1);

  Gens gs1;
  gs1.push_back(point());
  gs1.push_back(ray(A));
  gs1.push_back(ray(A - B));
  Poly ph1(2, Spec_Elem::EMPTY);
  ph1.add_gens(gs1);

  Gens gs2;
  gs2.push_back(point());
  gs2.push_back(ray(A));
  gs2.push_back(ray(A - 2*B));
  Poly ph2(2, Spec_Elem::EMPTY);
  ph2.add_gens(gs2);

  print_gens(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");

  ph2.widening_assign(ph1, Widen_Impl::BHRZ03);

  Poly known_result(2);
  known_result.add_con(A >= 0);
  known_result.add_con(B <= 0);

  bool ok = (ph2 == known_result);

  print_gens(ph2, "*** after ph2.widening_assign(ph1) ***");

  return ok;
}

bool
test16() {
  Var A(0);
  Var B(1);

  Gens gs1;
  gs1.push_back(point());
  gs1.push_back(ray(B));
  gs1.push_back(ray(A + 2*B));
  Poly ph1(2, Spec_Elem::EMPTY);
  ph1.add_gens(gs1);

  Gens gs2;
  gs2.push_back(point());
  gs2.push_back(ray(B));
  gs2.push_back(ray(A + B));
  Poly ph2(2, Spec_Elem::EMPTY);
  ph2.add_gens(gs2);

  print_gens(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");

  ph2.widening_assign(ph1, Widen_Impl::BHRZ03);

  Poly known_result(2);
  known_result.add_con(A >= 0);
  known_result.add_con(B >= 0);

  bool ok = (ph2 == known_result);

  print_gens(ph2, "*** after ph2.widening_assign(ph1) ***");

  return ok;
}

bool
test17() {
  Var A(0);
  Var B(1);

  Gens gs1;
  gs1.push_back(point());
  gs1.push_back(ray(B));
  gs1.push_back(ray(-A + 2*B));
  Poly ph1(2, Spec_Elem::EMPTY);
  ph1.add_gens(gs1);

  Gens gs2;
  gs2.push_back(point());
  gs2.push_back(ray(B));
  gs2.push_back(ray(-A + B));
  Poly ph2(2, Spec_Elem::EMPTY);
  ph2.add_gens(gs2);

  print_gens(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");

  ph2.widening_assign(ph1, Widen_Impl::BHRZ03);

  Poly known_result(2);
  known_result.add_con(A <= 0);
  known_result.add_con(B >= 0);

  bool ok = (ph2 == known_result);

  print_gens(ph2, "*** after ph2.widening_assign(ph1) ***");

  return ok;
}

bool
test18() {
  Var A(0);
  Var B(1);

  Gens gs1;
  gs1.push_back(point());
  gs1.push_back(ray(-B));
  gs1.push_back(ray(-A - 2*B));
  Poly ph1(2, Spec_Elem::EMPTY);
  ph1.add_gens(gs1);

  Gens gs2;
  gs2.push_back(point());
  gs2.push_back(ray(-B));
  gs2.push_back(ray(-A - B));
  Poly ph2(2, Spec_Elem::EMPTY);
  ph2.add_gens(gs2);

  print_gens(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");

  ph2.widening_assign(ph1, Widen_Impl::BHRZ03);

  Poly known_result(2);
  known_result.add_con(A <= 0);
  known_result.add_con(B <= 0);

  bool ok = (ph2 == known_result);

  print_gens(ph2, "*** after ph2.widening_assign(ph1) ***");

  return ok;
}

bool
test19() {
  Var A(0);
  Var B(1);

  Gens gs1;
  gs1.push_back(point());
  gs1.push_back(ray(-B));
  gs1.push_back(ray(A - 2*B));
  Poly ph1(2, Spec_Elem::EMPTY);
  ph1.add_gens(gs1);

  Gens gs2;
  gs2.push_back(point());
  gs2.push_back(ray(-B));
  gs2.push_back(ray(A - B));
  Poly ph2(2, Spec_Elem::EMPTY);
  ph2.add_gens(gs2);

  print_gens(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");

  ph2.widening_assign(ph1, Widen_Impl::BHRZ03);

  Poly known_result(2);
  known_result.add_con(A >= 0);
  known_result.add_con(B <= 0);

  bool ok = (ph2 == known_result);

  print_gens(ph2, "*** after ph2.widening_assign(ph1) ***");

  return ok;
}

bool
test20() {
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
  gs2.push_back(ray(-A + B));
  Poly ph2(2, Spec_Elem::EMPTY);
  ph2.add_gens(gs2);

  print_gens(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");

  ph2.widening_assign(ph1, Widen_Impl::BHRZ03);

  Poly known_result(2);
  known_result.add_con(B >= 0);

  bool ok = (ph2 == known_result);

  print_gens(ph2, "*** after ph2.widening_assign(ph1) ***");

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
  DO_TEST(test13);
  DO_TEST(test14);
  DO_TEST(test15);
  DO_TEST(test16);
  DO_TEST(test17);
  DO_TEST(test18);
  DO_TEST(test19);
  DO_TEST(test20);
END_MAIN
