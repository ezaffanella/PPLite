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

bool
test01() {
  Var A(0);

  Poly ph1(3);

  print_cons(ph1, "*** ph1 ***");

  ph1.expand_space_dim(A, 1);

  Poly kr(4);

  bool ok = (ph1 == kr);

  print_cons(ph1, "*** after ph1.expand_space_dim(A, 1) ***");

  return ok;
}

// Test with an empty polyhedron.
bool
test02() {
  //  Var A(0);
  Var B(1);

  Poly ph1(3, Spec_Elem::EMPTY);

  print_cons(ph1, "*** ph1 ***");

  ph1.expand_space_dim(B, 1);

  Poly kr(4, Spec_Elem::EMPTY);

  bool ok = (ph1 == kr);

  print_cons(ph1, "*** after ph1.expand_space_dim(B, 1) ***");

  return ok;
}

// Test trivial expansion.
bool
test03() {
  Var A(0);
  Var B(1);
  Var C(2);
  Var D(3);
  Var E(4);

  Poly ph1(2);
  ph1.add_con(A >= 0);
  ph1.add_con(A + B <= 2);

  print_cons(ph1, "*** ph1 ***");

  ph1.expand_space_dim(A, 0);

  Poly kr(2);
  kr.add_con(A >= 0);
  kr.add_con(A + B <= 2);

  bool ok = (ph1 == kr);

  print_cons(ph1, "*** after ph1.expand_space_dim(A, 0) ***");

  return ok;
}

// Test with given gens.
bool
test04() {
  Var A(0);
  Var B(1);
  Var C(2);

  Poly ph1(2, Spec_Elem::EMPTY);
  ph1.add_gen(point(A));
  ph1.add_gen(point(A + B));
  ph1.add_gen(point(B));

  print_gens(ph1, "*** ph1 ***");

  ph1.expand_space_dim(A, 1);

  Poly kr(3, Spec_Elem::EMPTY);
  kr.add_gen(point(A + C));
  kr.add_gen(point(A + B));
  kr.add_gen(point(A + B + C));
  kr.add_gen(point(B));
  kr.add_gen(point(B + C));

  bool ok = (ph1 == kr);

  print_gens(ph1, "*** after ph1.expand_space_dim(A, 1) ***");

  return ok;
}

// Test with given cons.
bool
test05() {
  Var A(0);
  Var B(1);
  Var C(2);

  Poly ph1(2);
  ph1.add_con(A >= 0);
  ph1.add_con(A + B <= 2);

  print_cons(ph1, "*** ph1 ***");

  ph1.expand_space_dim(A, 1);

  Poly kr(3);
  kr.add_con(A >= 0);
  kr.add_con(A + B <= 2);
  kr.add_con(C >= 0);
  kr.add_con(C + B <= 2);

  bool ok = (ph1 == kr);

  print_cons(ph1, "*** after ph1.expand_space_dim(A, 1) ***");

  return ok;
}

// Test using cons expanding 2 dimensions.
bool
test06() {
  Var A(0);
  Var B(1);
  Var C(2);
  Var D(3);

  Poly ph1(2);
  ph1.add_con(A >= 0);
  ph1.add_con(A + B <= 2);

  print_cons(ph1, "*** ph1 ***");

  ph1.expand_space_dim(A, 2);

  Poly kr(4);
  kr.add_con(A >= 0);
  kr.add_con(A + B <= 2);
  kr.add_con(C >= 0);
  kr.add_con(C + B <= 2);
  kr.add_con(D >= 0);
  kr.add_con(D + B <= 2);

  bool ok = (ph1 == kr);

  print_cons(ph1, "*** after ph1.expand_space_dim(A, 2) ***");

  return ok;
}

// Test using cons with equality con.
bool
test07() {
  Var A(0);
  Var B(1);
  Var C(2);
  Var D(3);
  Var E(4);

  Poly ph1(3);
  ph1.add_con(A <= 1);
  ph1.add_con(C == 1);
  ph1.add_con(A + B >= 1);
  ph1.add_con(B <= 1);

  print_cons(ph1, "*** ph1 ***");

  ph1.expand_space_dim(A, 1);
  ph1.expand_space_dim(C, 1);

  Poly kr(5);
  kr.add_con(A <= 1);
  kr.add_con(A + B >= 1);
  kr.add_con(C == 1);
  kr.add_con(E == 1);
  kr.add_con(B <= 1);
  kr.add_con(D <= 1);
  kr.add_con(D + B >= 1);

  bool ok = (ph1 == kr);

  print_cons(ph1,
                    "*** after ph1.expand_space_dim(A, 1);"
                    " ph1.expand_space_dim(C, 1) ***");

  return ok;
}

bool
test08() {
  Var A(0);
  Var B(1);
  Var C(2);

  Poly ph1(2, Spec_Elem::EMPTY);
  ph1.add_gen(point(A + 2*B));
  ph1.add_gen(point(A + 3*B));
  ph1.add_gen(point(A + 4*B));

  print_gens(ph1, "*** ph1 ***");

  ph1.expand_space_dim(B, 1);

  Poly kr(3, Spec_Elem::EMPTY);
  kr.add_gen(point(A + 2*B + 2*C));
  kr.add_gen(point(A + 2*B + 3*C));
  kr.add_gen(point(A + 2*B + 4*C));
  kr.add_gen(point(A + 3*B + 2*C));
  kr.add_gen(point(A + 3*B + 3*C));
  kr.add_gen(point(A + 3*B + 4*C));
  kr.add_gen(point(A + 4*B + 2*C));
  kr.add_gen(point(A + 4*B + 3*C));
  kr.add_gen(point(A + 4*B + 4*C));

  bool ok = (ph1 == kr);

  print_gens(ph1, "*** after ph1.expand_space_dim(A, 2) ***");

  return ok;
}

bool
test09() {
  Var A(0);
  Var B(1);
  Var C(2);
  Var D(3);

  Poly ph1(2, Topol::NNC);
  ph1.add_con(A - B > 2);
  ph1.add_con(A + 2*B < 6);
  ph1.add_con(B < 6);

  print_cons(ph1, "*** ph1 ***");

  ph1.expand_space_dim(B, 2);

  Poly kr(4, Topol::NNC);
  kr.add_con(A - B > 2);
  kr.add_con(A + 2*B < 6);
  kr.add_con(B < 6);
  kr.add_con(A - C > 2);
  kr.add_con(A + 2*C < 6);
  kr.add_con(C < 6);
  kr.add_con(A - D > 2);
  kr.add_con(A + 2*D < 6);
  kr.add_con(D < 6);

  bool ok = (ph1 == kr);

  print_cons(ph1, "*** after ph1.expand_space_dim(B, 2) ***");

  return ok;
}

bool
test10() {
  Var A(0);
  Var B(1);
  Var C(2);
  Var D(3);

  Poly ph1(2, Topol::NNC, Spec_Elem::EMPTY);
  ph1.add_gen(point(A));
  ph1.add_gen(closure_point(A + B));
  ph1.add_gen(ray(A - B));

  print_gens(ph1, "*** ph1 ***");

  ph1.expand_space_dim(A, 2);

  Poly kr(4, Topol::NNC, Spec_Elem::EMPTY);
  kr.add_gen(point(A + C + D));
  kr.add_gen(ray(A -B + C + D));
  kr.add_gen(closure_point(A + C + 2*D));
  kr.add_gen(closure_point(A + 2*C + D));
  kr.add_gen(closure_point(A + 2*C + 2*D));
  kr.add_gen(closure_point(A + B + C + D));
  kr.add_gen(closure_point(2*A + C + D));
  kr.add_gen(closure_point(2*A + C + 2*D));
  kr.add_gen(closure_point(2*A + 2*C + D));

  bool ok = (ph1 == kr);

  print_gens(ph1, "*** after ph1.expand_space_dim(A, 2) ***");

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
