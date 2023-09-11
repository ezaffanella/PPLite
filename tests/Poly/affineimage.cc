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
  Var C(2);

  Poly ph(3);
  ph.add_con(C == -2);
  ph.add_con(A == 0);

  print_cons(ph, "*** ph ***");

  ph.affine_image(B, A, 2, 1);

  Poly known_result(3, Spec_Elem::EMPTY);
  known_result.add_gen(point(2*B - 2*C));

  bool ok = (ph == known_result);

  print_gens(ph, "*** ph after ph.affine_image(B, A, 2, 1) ***");

  return ok;
}

bool
test02() {
  Var A(0);
  Var B(1);

  Poly ph(2);
  ph.add_con(A >= B);
  ph.add_con(B >= 0);
  ph.add_con(A <= 3);

  print_cons(ph, "*** ph ***");

  ph.affine_image(A, A+B, 1);

  Poly known_result(2);
  known_result.add_con(A - 2*B >= 1);
  known_result.add_con(B >= 0);
  known_result.add_con(A - B <= 4);

  bool ok = (ph == known_result);

  print_cons(ph, "*** ph after ph.affine_image(A, A+B, 1) ***");

  return ok;
}

bool
test03() {
  Var A(0);
  Var B(1);

  Gens gs;
  gs.push_back(point());
  gs.push_back(ray(A));
  gs.push_back(ray(B));
  Poly ph(2, Spec_Elem::EMPTY);
  ph.add_gens(gs);

  print_gens(ph, "*** ph ***");

  ph.affine_image(A, A, 1, 2);

  Poly known_result(2);
  known_result.add_con(2*A >= 1);
  known_result.add_con(B >= 0);

  bool ok = (ph == known_result);

  print_gens(ph, "*** after ph.affine_image(A, A, 1, 2) ***");

  return ok;
}

bool
test04() {
  Var A(0);
  Var B(1);

  Poly ph(2, Spec_Elem::EMPTY);
  ph.add_gen(point(A));

  print_cons(ph, "*** ph ***");

  ph.affine_image(A, B, 2, -3);

  Poly known_result(2, Spec_Elem::EMPTY);
  known_result.add_gen(point(-2*A, 3));

  bool ok = (ph == known_result);

  print_gens(ph, "*** ph after ph.affine_image(A, B, 2, -3) ***");

  return ok;
}

bool
test05() {
  Var A(0);
  Var B(1);

  Poly ph(2);
  ph.add_con(A >= 2);
  ph.add_con(A <= 3);
  ph.add_con(B >= 1);
  ph.add_con(2*A >= B);

  print_cons(ph, "*** ph ***");

  ph.affine_image(B, A-B, 2, -3);

  Poly known_result(2, Spec_Elem::EMPTY);
  known_result.add_gen(point(2*A));
  known_result.add_gen(point(2*A - B));
  known_result.add_gen(point(9*A + B, 3));
  known_result.add_gen(point(9*A - 4*B, 3));

  bool ok = (ph == known_result);

  print_gens(ph, "*** ph after ph.affine_image(B, A-B, 2, -3) ***");

  return ok;
}

bool
test06() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::EMPTY);

  print_cons(ph1, "*** ph1 ***");

  ph1.affine_image(A, 2*A + B, 1);

  Poly known_result(2, Spec_Elem::EMPTY);

  bool ok = (ph1 == known_result);

  print_cons(ph1, "*** after ph1.affine_image(A, 2*A + B, 1) ***");

  return ok;
}

bool
test07() {
  Var A(0);
  Var B(1);

  Gens gs;
  gs.push_back(point());
  gs.push_back(point(A));
  gs.push_back(point(B));
  gs.push_back(point(A + B));

  Poly ph(2, Spec_Elem::EMPTY);
  ph.add_gens(gs);

  print_gens(ph, "*** ph ***");

  ph.affine_image(A, -A, -1, -1);

  Gens known_gs;
  known_gs.push_back(point(A));
  known_gs.push_back(point(2*A));
  known_gs.push_back(point(A + B));
  known_gs.push_back(point(2*A + B));
  Poly known_result(2, Spec_Elem::EMPTY);
  known_result.add_gens(known_gs);

  bool ok = (ph == known_result);

  print_gens(ph, "*** after ph.affine_image(A, -A, -1, -1) ***");

  return ok;
}

bool
test08() {
  Var A(0);
  Var B(1);

  Poly ph(2);
  ph.add_con(A >= 0);
  ph.add_con(B >= 0);
  Poly copy_ph(ph);

  print_cons(ph, "*** ph ***");

  ph.affine_image(A, A, 1);
  copy_ph.affine_image(A, -A, -1, -1);

  bool ok = (ph == copy_ph);

  print_gens(ph, "*** after ph.affine_image(A, A, 1) ***");
  print_gens(copy_ph, "*** after copy_ph.affine_image(A, -A, -1, -1) ***");

  return ok;
}

bool
test09() {
  Var A(0);
  Var B(1);

  Poly ph(2);
  ph.add_con(A >= 0);
  ph.add_con(B >= 0);
  Poly copy_ph(ph);

  print_cons(ph, "*** ph ***");

  ph.affine_image(B, A, 1);
  copy_ph.affine_image(B, -A, -1, -1);

  bool ok = (ph == copy_ph);

  print_gens(ph, "*** after ph.affine_image(B, A, 1) ***");
  print_gens(copy_ph, "*** after copy_ph.affine_image(B, -A, -1, -1) ***");

  return ok;
}

bool
test10() {
  Var A(0);
  Var B(1);

  Poly ph(2, Spec_Elem::EMPTY);
  ph.add_gen(point(A));
  print_gens(ph, "*** ph ***");

  ph.affine_image(A, Linear_Expr());
  print_gens(ph, "*** after ph.affine_image(A, 0) ***");

  Poly known_value(2, Spec_Elem::EMPTY);
  known_value.add_gen(point());
  bool ok = (ph == known_value);

  return ok;
}

bool
test11() {
  Var A(0);
  Var B(1);
  Var C(2);

  Poly ph(3);
  print_cons(ph, "*** ph ***");

  ph.affine_image(A, C);
  print_cons(ph, "*** after ph.affine_image(A, C) ***");

  Poly known_value(3);
  known_value.add_con(A == C);
  bool ok = (ph == known_value);

  return ok;
}


bool
test12() {
  Var A(0);
  Var B(1);
  Var C(2);

  Poly ph(3);
  ph.add_con(A - B == 1);
  print_cons(ph, "*** ph ***");

  ph.affine_image(A, C);
  print_cons(ph, "*** after ph.affine_image(A, C) ***");

  Poly known_value(3);
  known_value.add_con(A == C);
  bool ok = (ph == known_value);

  return ok;
}

bool
test13() {
  Var A(0);
  Var B(1);

  Poly p1(2, Topol::NNC);
  p1.add_con(B == 0);
  p1.add_con(A >= 0);

  Poly p2(2, Topol::NNC);
  p2.add_con(B == 0);
  p2.add_con(-A > 0);

  print_cons(p1, "*** p1 ***");
  print_cons(p2, "*** p2 ***");

  p1.affine_image(B, A, 2);
  p2.affine_image(B, A);

  print_cons(p1, "*** p1.affine_image(B, A, 2) ***");
  print_cons(p2, "*** p2.affine_image(B, A) ***");

  p1.poly_hull_assign(p2);

  Poly known_result(2, Topol::NNC);
  known_result.add_con(A - B >= -2);
  known_result.add_con(-A + B >= 0);

  bool ok = (p1 == known_result);

  print_cons(p1, "*** p1.upper_bound_assign(p2) ***");

  return ok;
}

bool
test14_aux(dim_type dim) {
  Poly ph(dim);
  for (auto i = 0; i < dim; ++i) {
    ph.add_con(Var(i) >= 0);
    ph.add_con(Var(i) <= 2);
  }
  ph.minimize();

  {
    Clock clock;
    ph.affine_image(Var(0), Var(1), 10);
    ph.minimize();
    nout << "Time for dim = " << dim << " : ";
    clock.print_elapsed(nout);
    nout << "\n";
  }

  Poly kr(dim);
  kr.add_con(Var(0) >= 10);
  kr.add_con(Var(0) - Var(1) == 10);
  for (auto i = 1; i < dim; ++i) {
    kr.add_con(Var(i) >= 0);
    kr.add_con(Var(i) <= 2);
  }
  return (ph == kr);
}

bool
test14() {
  const dim_type min_dim = 8;
  const dim_type max_dim = 12;
  for (auto i = min_dim; i <= max_dim; ++i) {
    if (!test14_aux(i))
      return false;
  }
  return true;
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
END_MAIN
