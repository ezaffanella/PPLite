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
  Var A(0);
  Var B(1);

  Poly ph1(3);
  print_gens(ph1, "*** ph1 ***");

  Vars_Set to_fold;
  to_fold.insert(A);
  ph1.fold_space_dims(to_fold, B);

  Poly known_result(2);

  bool ok = (ph1 == known_result);

  print_gens(ph1, "*** after folding {A} into B ***");

  return ok;
}

bool
test02() {
  Var A(0);
  Var B(1);

  Poly ph1(3, Spec_Elem::EMPTY);
  print_cons(ph1, "*** ph1 ***");

  Vars_Set to_fold;
  to_fold.insert(A);

  ph1.fold_space_dims(to_fold, B);

  Poly known_result(2, Spec_Elem::EMPTY);

  bool ok = (ph1 == known_result);

  print_cons(ph1, "*** after folding {A} into B ***");

  return ok;
}

bool
test03() {
  Var A(0);
  Var B(1);
  Var C(2);

  Poly ph1(3);
  ph1.add_con(A >= 0);
  ph1.add_con(A + B + C <= 2);

  print_cons(ph1, "*** ph1 ***");

  Vars_Set to_fold;

  ph1.fold_space_dims(to_fold, B);

  Poly known_result(3);
  known_result.add_con(A >= 0);
  known_result.add_con(A + B + C <= 2);

  bool ok = (ph1 == known_result);

  print_cons(ph1, "*** after folding {} into B ***");

  return ok;
}

bool
test04() {
  Var A(0);
  Var B(1);

  Poly ph1(2);
  ph1.add_con(A >= 1);
  ph1.add_con(A <= 3);
  ph1.add_con(B >= 7);
  ph1.add_con(B <= 12);

  print_cons(ph1, "*** ph1 ***");

  Vars_Set to_fold;
  to_fold.insert(A);

  ph1.fold_space_dims(to_fold, B);

  Poly known_result(1);
  known_result.add_con(A >= 1);
  known_result.add_con(A <= 12);

  bool ok = (ph1 == known_result);

  print_cons(ph1, "*** after folding {A} into B ***");

  return ok;
}

bool
test05() {
  Var A(0);
  Var B(1);
  Var C(2);

  Poly ph1(3, Spec_Elem::EMPTY);
  ph1.add_gen(point(A + 2*B + 2*C));
  ph1.add_gen(point(A + 2*B + 3*C));
  ph1.add_gen(point(A + 2*B + 4*C));
  ph1.add_gen(point(A + 3*B + 2*C));
  ph1.add_gen(point(A + 3*B + 3*C));
  ph1.add_gen(point(A + 3*B + 4*C));
  ph1.add_gen(point(A + 4*B + 2*C));
  ph1.add_gen(point(A + 4*B + 3*C));
  ph1.add_gen(point(A + 4*B + 4*C));

  print_gens(ph1, "*** ph1 ***");

  Vars_Set to_fold;
  to_fold.insert(C);

  ph1.fold_space_dims(to_fold, B);

  Poly known_result(2, Spec_Elem::EMPTY);
  known_result.add_gen(point(A + 2*B));
  known_result.add_gen(point(A + 3*B));
  known_result.add_gen(point(A + 4*B));

  bool ok = (ph1 == known_result);

  print_gens(ph1, "*** after folding {C} into B ***");

  return ok;
}

bool
test06() {
  Var A(0);
  Var B(1);
  Var C(2);

  Poly ph1(3);
  ph1.add_con(A >= 1);
  ph1.add_con(A <= 3);
  ph1.add_con(B >= 7);
  ph1.add_con(B <= 12);
  ph1.add_con(C == 15);

  print_cons(ph1, "*** ph1 ***");

  Vars_Set to_fold;
  to_fold.insert(A);
  to_fold.insert(B);

  ph1.fold_space_dims(to_fold, C);

  Poly known_result(1);
  known_result.add_con(A >= 1);
  known_result.add_con(A <= 15);

  bool ok = (ph1 == known_result);

  print_cons(ph1, "*** after folding {A,B} into C ***");

  return ok;
}

bool
test07() {
  Var A(0);
  Var B(1);
  Var C(2);

  Poly ph1(3, Spec_Elem::EMPTY);
  ph1.add_gen(point(A));
  ph1.add_gen(ray(A + B));
  ph1.add_gen(ray(A + 2*C));

  print_gens(ph1, "*** ph1 ***");

  Vars_Set to_fold;
  to_fold.insert(C);

  ph1.fold_space_dims(to_fold, B);

  Poly known_result(2, Spec_Elem::EMPTY);
  known_result.add_gen(point(A));
  known_result.add_gen(ray(A));
  known_result.add_gen(ray(A + B));
  known_result.add_gen(ray(A + 2*B));

  bool ok = (ph1 == known_result);

  print_gens(ph1, "*** after folding {C} into B ***");

  return ok;
}

bool
test08() {
  Var A(0);
  Var B(1);
  Var C(2);
  Var D(3);

  Poly ph1(4);
  ph1.add_con(A >= 0);
  ph1.add_con(A + B <= 2);
  ph1.add_con(C >= 0);
  ph1.add_con(C + B <= 2);
  ph1.add_con(D >= 0);
  ph1.add_con(D + B <= 2);

  print_cons(ph1, "*** ph1 ***");

  Vars_Set to_fold;
  to_fold.insert(C);
  to_fold.insert(D);

  ph1.fold_space_dims(to_fold, A);

  Poly known_result(2);
  known_result.add_con(A >= 0);
  known_result.add_con(A + B <= 2);

  bool ok = (ph1 == known_result);

  print_cons(ph1, "*** after folding {C,D} into A ***");

  return ok;
}

bool
test09() {
  Var A(0);
  Var B(1);
  Var C(2);
  Var D(3);

  Poly ph1(4);
  ph1.add_con(A >= 0);
  ph1.add_con(B >= 0);
  ph1.add_con(A + B <= 2);
  ph1.add_con(C >= 0);
  ph1.add_con(C + B <= 2);
  ph1.add_con(D >= 0);
  ph1.add_con(D + B <= 2);

  print_cons(ph1, "*** ph1 ***");

  Vars_Set to_fold;
  to_fold.insert(B);
  to_fold.insert(D);

  ph1.fold_space_dims(to_fold, C);

  Poly known_result(2);
  known_result.add_con(A >= 0);
  known_result.add_con(A <= 2);
  known_result.add_con(B >= 0);
  known_result.add_con(B <= 2);

  bool ok = (ph1 == known_result);

  print_cons(ph1, "*** after folding {B,D} into C ***");

  return ok;
}

bool
test10() {
  Var A(0);
  Var B(1);

  Poly ph1(3, Topol::NNC, Spec_Elem::EMPTY);

  print_cons(ph1, "*** ph1 ***");

  Vars_Set to_fold;
  to_fold.insert(A);

  ph1.fold_space_dims(to_fold, B);

  Poly known_result(2, Topol::NNC, Spec_Elem::EMPTY);

  bool ok = (ph1 == known_result);

  print_cons(ph1, "*** after folding {A} into B ***");

  return ok;
}

bool
test11() {
  Var A(0);
  Var B(1);
  Var C(2);

  Poly ph1(3, Topol::NNC);
  ph1.add_con(A >= 0);
  ph1.add_con(A + B + C < 2);

  print_cons(ph1, "*** ph1 ***");

  Vars_Set to_fold;

  ph1.fold_space_dims(to_fold, B);

  Poly known_result(3, Topol::NNC);
  known_result.add_con(A >= 0);
  known_result.add_con(A + B + C < 2);

  bool ok = (ph1 == known_result);

  print_cons(ph1, "*** after folding {} into B ***");

  return ok;
}

bool
test12() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Topol::NNC);
  ph1.add_con(A > 1);
  ph1.add_con(A < 3);
  ph1.add_con(B > 7);
  ph1.add_con(B < 12);

  print_cons(ph1, "*** ph1 ***");

  Vars_Set to_fold;
  to_fold.insert(A);

  ph1.fold_space_dims(to_fold, B);

  Poly known_result(1, Topol::NNC);
  known_result.add_con(A > 1);
  known_result.add_con(A < 12);

  bool ok = (ph1 == known_result);

  print_cons(ph1, "*** after folding {A} into B ***");

  return ok;
}

bool
test13() {
  Var A(0);
  Var B(1);
  Var C(2);

  Poly ph1(3, Topol::NNC);
  ph1.add_con(A > 1);
  ph1.add_con(A <= 3);
  ph1.add_con(B > 7);
  ph1.add_con(B < 12);
  ph1.add_con(C == 15);

  print_cons(ph1, "*** ph1 ***");
  print_gens(ph1, "*** ph1 ***");

  Vars_Set to_fold;
  to_fold.insert(A);
  to_fold.insert(B);

  ph1.fold_space_dims(to_fold, C);

  Poly known_result(1, Topol::NNC);
  known_result.add_con(A > 1);
  known_result.add_con(A <= 15);

  bool ok = (ph1 == known_result);

  print_cons(ph1, "*** after folding {A,B} into C ***");

  return ok;
}

bool
test14() {
  Var A(0);
  Var B(1);
  Var C(2);
  Var D(3);

  Poly ph1(4, Topol::NNC);
  ph1.add_con(A > 0);
  ph1.add_con(A + B < 2);
  ph1.add_con(C > 0);
  ph1.add_con(C + B < 2);
  ph1.add_con(D > 0);
  ph1.add_con(D + B <= 2);

  print_cons(ph1, "*** ph1 ***");

  Vars_Set to_fold;
  to_fold.insert(C);
  to_fold.insert(D);

  ph1.fold_space_dims(to_fold, A);

  Poly known_result(2, Topol::NNC);
  known_result.add_con(A > 0);
  known_result.add_con(A + B <= 2);

  bool ok = (ph1 == known_result);

  print_cons(ph1, "*** after folding {C,D} into A ***");

  return ok;
}

bool
test15() {
  Var A(0);
  Var B(1);
  Var C(2);
  Var D(3);

  Poly ph1(4, Topol::NNC);
  ph1.add_con(A >= 0);
  ph1.add_con(B > 0);
  ph1.add_con(A + B <= 2);
  ph1.add_con(C > 0);
  ph1.add_con(C + B <= 2);
  ph1.add_con(D > 0);
  ph1.add_con(D + B <= 2);

  print_cons(ph1, "*** ph1 ***");

  Vars_Set to_fold;
  to_fold.insert(B);
  to_fold.insert(D);

  ph1.fold_space_dims(to_fold, C);

  Poly known_result(2, Topol::NNC);
  known_result.add_con(A >= 0);
  known_result.add_con(A < 2);
  known_result.add_con(B > 0);
  known_result.add_con(B < 2);

  bool ok = (ph1 == known_result);

  print_cons(ph1, "*** after folding {B,D} into C ***");

  return ok;
}

bool
test16() {
  Var A(0);
  Var B(1);

  Poly ph(2, Topol::NNC);
  ph.add_con(A <= 2);
  ph.add_con(B > 0);
  ph.add_con(B < 1);

  print_cons(ph, "*** ph ***");
  print_gens(ph, "*** ph ***");

  Vars_Set to_fold;
  to_fold.insert(B);

  ph.fold_space_dims(to_fold, A);

  Poly known_result(1, Topol::NNC);
  known_result.add_con(A <= 2);

  bool ok = (ph == known_result);

  print_cons(ph, "*** after folding {B} into A ***");

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
END_MAIN
