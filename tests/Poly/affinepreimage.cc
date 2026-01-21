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

  Gens gs;
  gs.push_back(point());
  gs.push_back(ray(A + B));
  gs.push_back(point(A));

  Poly ph(2, Spec_Elem::EMPTY);
  ph.add_gens(gs);
  print_gens(ph, "*** ph ***");
  ph.affine_preimage(A, A, 2);
  print_gens(ph, "*** ph after ph.affine_preimage(A, A, 2) ***");

  Gens kr_gs;
  kr_gs.push_back(point(-2*A));
  kr_gs.push_back(ray(A + B));
  kr_gs.push_back(point(-A));
  Poly kr(2, Spec_Elem::EMPTY);
  kr.add_gens(kr_gs);
  return (ph == kr);
}

bool
test02() {
  Var A(0);
  Var B(1);
  Var C(2);

  Poly ph(3);
  ph.add_con(A + C == 0);
  ph.add_con(A + B >= 0);
  ph.add_con(A + B >= 1);

  print_cons(ph, "*** ph ***");

  ph.affine_preimage(A, A + B);

  Poly known_result(3);
  known_result.add_con(A + B + C == 0);
  known_result.add_con(A + 2*B >= 0);
  known_result.add_con(A + 2*B >= 1);

  bool ok = (ph == known_result);

  print_cons(ph, "*** ph after ph.affine_preimage(A, A+B) ***");

  return ok;
}

bool
test03() {
  Var A(0);
  Var B(1);

  Poly ph(2);
  ph.add_con(A >= 0);
  ph.add_con(B >= 0);
  ph.add_con(A + B >= 3);

  print_cons(ph, "*** ph ***");

  ph.affine_preimage(A, B, 1);

  Poly known_result(2, Spec_Elem::EMPTY);
  known_result.add_gen(point(B));
  known_result.add_gen(line(A));
  known_result.add_gen(ray(B));

  bool ok = (ph == known_result);

  print_gens(ph, "*** ph after ph.affine_preimage(A, B, 1) ***");

  return ok;
}

bool
test04() {
  Var x(0);
  Var y(1);

  Gens gs;
  gs.push_back(point(x + y));
  gs.push_back(ray(x + 2*y));
  gs.push_back(ray(x));

  Poly ph(2, Spec_Elem::EMPTY);
  ph.add_gens(gs);
  print_gens(ph, "*** ph ***");
  ph.affine_preimage(x, y, 1);

  Poly known_result(2);
  known_result.add_con(y >= 1);

  bool ok = (ph == known_result);

  print_cons(ph, "*** after affine_preimage(x, y, 1) ***");

  return ok;
}

bool
test05() {
  Var A(0);
  Var B(1);

  Poly ph(2);
  ph.add_con(A >= 2);
  ph.add_con(B >= 0);

  print_cons(ph, "*** ph ***");

  ph.affine_preimage(A, A, 1, 2);

  Poly known_result(2);
  known_result.add_con(A >= 3);
  known_result.add_con(B >= 0);

  bool ok = (ph == known_result);

  print_cons(ph, "*** after ph.affine_preimage(A, A, 1, 2) ***");

  return ok;
}

bool
test06() {
  Var A(0);
  Var B(1);
  Poly ph(2);
  ph.add_con(A >= 0);
  ph.add_con(B >= 2);

  print_cons(ph, "*** ph ***");

  ph.affine_preimage(B, A, 1, 2);

  Poly known_result(2);
  known_result.add_con(A >= 3);

  bool ok = (ph == known_result);

  print_cons(ph, "*** after ph.affine_preimage(B, A, 1, 2) ***");

  return ok;
}

bool
test07() {
  Var A(0);
  Var B(1);

  Gens gs;
  gs.push_back(point());
  gs.push_back(ray(A));
  gs.push_back(ray(A + B));
  Poly ph(2, Spec_Elem::EMPTY);
  ph.add_gens(gs);
  print_gens(ph, "*** ph ***");
  ph.affine_preimage(B, A - B, 0, -1);

  Poly known_result(2);
  known_result.add_con(A - B <= 0);
  known_result.add_con(2*A - B >= 0);

  bool ok = (ph == known_result);

  print_cons(ph, "*** after ph.affine_preimage(B, A - B, 0, -1) ***");

  return ok;
}

bool
test08() {
  Var A(0);
  Var B(1);

  Poly ph(2, Spec_Elem::EMPTY);

  print_cons(ph, "*** ph ***");

  ph.affine_preimage(A, 2*A + B, 1);

  Poly known_result(2, Spec_Elem::EMPTY);

  bool ok = (ph == known_result);

  print_cons(ph, "*** after ph.affine_preimage(A, 2*A + B, 1) ***");

  return ok;
}

bool
test09() {
  Var A(0);
  Var B(1);

  Poly ph(2);
  ph.add_con(A >= 0);
  ph.add_con(B >= 0);
  ph.add_con(A <= 1);
  ph.add_con(B <= 1);

  ph.affine_preimage(A, -A, -1, -1);

  Poly known_result(2);
  known_result.add_con(A <= 0);
  known_result.add_con(B <= 1);
  known_result.add_con(A >= -1);
  known_result.add_con(B >= 0);

  bool ok = (ph == known_result);

  print_cons(ph, "*** after ph.affine_preimage(A, -A, -1, -1) ***");

  return ok;
}

bool
test10() {
  Var A(0);
  Var B(1);

  Poly ph(2);
  ph.add_con(A >= 0);
  ph.add_con(B >= 0);
  ph.add_con(A <= 1);
  ph.add_con(B <= 1);

  print_cons(ph, "*** ph ***");

  ph.affine_preimage(B, -A, -1, -1);

  Poly known_result(2);
  known_result.add_con(A == 0);

  bool ok = (ph == known_result);

  print_cons(ph, "*** after ph.affine_preimage(B, -A, -1, -1) ***");

  return ok;
}

bool
test11() {
  Var A(0);
  Var B(1);

  Poly ph(2);
  ph.add_con(A >= 0);
  ph.add_con(B >= 0);
  Poly copy_ph(ph);

  print_cons(ph, "*** ph ***");

  ph.affine_preimage(A, A, 1);
  copy_ph.affine_preimage(A, -A, -1, -1);

  bool ok = (ph == copy_ph);

  print_gens(ph, "*** after ph.affine_preimage(A, A, 1) ***");
  print_gens(copy_ph, "*** after copy_ph.affine_preimage(A, -A, -1, -1) ***");
  return ok;
}

bool
test12() {
  Var A(0);
  Var B(1);

  Poly ph(2);
  ph.add_con(A >= 0);
  ph.add_con(B >= 0);
  Poly copy_ph(ph);

  print_cons(ph, "*** ph ***");

  ph.affine_preimage(B, A, 1);
  copy_ph.affine_preimage(B, -A, -1, -1);

  bool ok = (ph == copy_ph);

  print_gens(ph, "*** after ph.affine_preimage(B, A, 1) ***");
  print_gens(copy_ph, "*** after copy_ph.affine_preimage(B, -A, -1, -1) ***");
  return ok;
}

bool
test13() {
  Var A(0);
  Var B(1);

  Gens gs;
  gs.push_back(point());
  gs.push_back(ray(A + B));
  gs.push_back(point(A));
  Poly ph(2, Spec_Elem::EMPTY);
  ph.add_gens(gs);
  print_gens(ph, "*** ph ***");

  ph.affine_preimage(B, Linear_Expr(), 1);
  print_gens(ph, "*** ph after ph.affine_preimage(B, Linear_Expr(), 1) ***");

  Gens kr_gs;
  kr_gs.push_back(point(A));
  kr_gs.push_back(point(2*A));
  kr_gs.push_back(line(B));
  Poly kr(2, Spec_Elem::EMPTY);
  kr.add_gens(kr_gs);

  return (ph == kr);
}

bool
test14() {
  Var A(0);
  Var B(1);

  Gens gs;
  gs.push_back(point());
  gs.push_back(point(0*A + 3*B));
  gs.push_back(point(3*A + 0*B));
  gs.push_back(point(3*A + 3*B));

  Poly ph(2, Spec_Elem::EMPTY);
  ph.add_gens(gs);

  Poly known_result = ph;

  print_gens(ph, "*** ph before ph.affine_image(A, A + 2*B, 4) ***");

  ph.affine_image(A, A + 2*B, 4);

  print_gens(ph, "*** ph after ph.affine_image(A, A + 2*B, 4) ***");

  ph.affine_preimage(A, A + 2*B, 4);

  bool ok = (ph == known_result);

  print_gens(ph, "*** ph after ph.affine_preimage(A, A + 2*B, 4) ***");

  return ok;
}

bool
test15() {
  Var A(0);
  Var B(1);

  Gens gs;
  gs.push_back(point());
  gs.push_back(point(0*A + 3*B));
  gs.push_back(point(3*A + 0*B));
  gs.push_back(point(3*A + 3*B));
  Poly ph(2, Spec_Elem::EMPTY);
  ph.add_gens(gs);

  print_gens(ph, "*** ph before ph.affine_image(A, B) ***");
  ph.affine_image(A, B);
  print_gens(ph, "*** ph after ph.affine_image(A, B) ***");
  ph.affine_preimage(A, B);
  print_gens(ph, "*** ph after ph.affine_preimage(A, B) ***");

  Poly known_result(2);
  known_result.add_con(B >= 0);
  known_result.add_con(B <= 3);

  return (ph == known_result);
}

bool
test16_aux(dim_type dim) {
  Poly ph(dim);
  for (auto i = 0; i < dim; ++i) {
    ph.add_con(Var(i) >= 0);
    ph.add_con(Var(i) <= 2);
  }
  ph.minimize();

  {
    Clock clock;
    ph.affine_preimage(Var(0), Var(1), 1);
    ph.minimize();
    nout << "Time for dim = " << dim << " : ";
    clock.print_elapsed(nout);
    nout << "\n";
  }

  Poly kr(dim);
  kr.add_con(Var(1) >= 0);
  kr.add_con(Var(1) <= 1);
  for (auto i = 2; i < dim; ++i) {
    kr.add_con(Var(i) >= 0);
    kr.add_con(Var(i) <= 2);
  }
  return (ph == kr);
}

bool
test16() {
  const dim_type min_dim = 8;
  const dim_type max_dim = 12;
  for (auto i = min_dim; i <= max_dim; ++i) {
    if (!test16_aux(i))
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
  DO_TEST(test15);
  DO_TEST(test16);
END_MAIN
