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
  Var z(2);
  Var w(3);
  Cons cs1 { x >= 0, y >= 0, x - y >= 0 };
  Poly ph(2);
  ph.add_cons(cs1);

  print_cons(ph, "*** ph before ***");

  Cons cs2 { x >= 1, y >= 1, x - y >= -1 };
  Poly qh(2);
  qh.add_cons(cs2);

  Poly copy_ph = ph;

  ph.concatenate_assign(qh);

  copy_ph.add_space_dims(2);
  copy_ph.add_con(z >= 1);
  copy_ph.add_con(w >= 1);
  copy_ph.add_con(z - w >= -1);

  bool ok = (ph == copy_ph);

  print_cons(ph, "*** concatenate_assign ***");
  print_cons(copy_ph, "*** embed + renaming + add_cons ***");

  return ok;
}

bool
test02() {
  Var x(0);
  Var y(1);

  Poly ph(2, Spec_Elem::EMPTY);

  print_cons(ph, "*** ph ***");

  Cons cs { x >= y, x >= 2 };
  Poly qh(2);
  qh.add_cons(cs);

  print_cons(qh, "*** qh ***");

  ph.concatenate_assign(qh);

  Poly known_result(4, Spec_Elem::EMPTY);

  bool ok = (ph == known_result);

  print_cons(ph, "*** after concatenate_assign(qh) ***");

  return ok;
}

bool
test03() {
  Var x(0);
  Var y(1);

  Poly ph;

  print_cons(ph, "*** ph ***");

  Cons cs { x - 3 >= y, y >= 0 };
  Poly qh(2);
  qh.add_cons(cs);

  print_cons(qh, "*** qh ***");

  ph.concatenate_assign(qh);

  bool ok = (ph == qh);

  print_cons(ph, "*** after concatenate_assign(qh) ***");

  return ok;
}

bool
test04() {
  Var A(0);
  Var B(1);
  Var C(2);
  Var D(3);

  Gens gs1;
  gs1.push_back(point());
  gs1.push_back(point(2*A));
  gs1.push_back(point(2*B));
  gs1.push_back(point(2*A + 2*B));
  Poly ph1(2, Spec_Elem::EMPTY);
  ph1.add_gens(gs1);

  Gens gs2;
  gs2.push_back(point(2*A));
  gs2.push_back(point(2*A + 3*B));
  Poly ph2(2, Spec_Elem::EMPTY);
  ph2.add_gens(gs2);

  print_gens(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");

  ph1.concatenate_assign(ph2);

  Poly known_result(4);
  known_result.add_con(A >= 0);
  known_result.add_con(B >= 0);
  known_result.add_con(A <= 2);
  known_result.add_con(B <= 2);
  known_result.add_con(C == 2);
  known_result.add_con(D >= 0);
  known_result.add_con(D <= 3);

  bool ok = (ph1 == known_result);

  print_cons(ph1, "*** after concatenate_assign(ph2) ***");

  return ok;
}

bool
test05() {
  Var A(0);
  Var B(1);

  Poly ph1(2);
  ph1.add_con(B >= 0);
  ph1.add_con(A - B >= 0);

  Poly ph2;

  Poly known_result = ph1;

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  ph1.concatenate_assign(ph2);

  bool ok = (ph1 == known_result);

  print_cons(ph1, "*** after ph1.concatenate_assign(ph2) ***");

  return ok;
}

bool
test06() {
  Var A(0);
  Var B(1);

  Poly ph1(1);
  ph1.gens();
  ph1.add_con(A >= 0);

  Poly ph2(1);
  ph2.gens();
  ph2.add_con(A == 2);

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  ph1.concatenate_assign(ph2);

  Poly known_result(2);
  known_result.add_con(A >= 0);
  known_result.add_con(B == 2);

  bool ok = (ph1 == known_result);

  print_cons(ph1, "*** after ph1.concatenate_assign(ph2) ***");

  return ok;
}

bool
test07() {
  Var A(0);
  Var B(1);

  Poly ph1(1, Spec_Elem::EMPTY);
  ph1.add_gen(point());
  ph1.minimize();
  ph1.add_gen(ray(A));

  Poly ph2(1, Spec_Elem::EMPTY);
  ph2.add_gen(point(2*A));
  ph2.minimize();

  print_gens(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");

  ph1.concatenate_assign(ph2);

  Poly known_result(2);
  known_result.add_con(A >= 0);
  known_result.add_con(B == 2);

  bool ok = (ph1 == known_result);

  print_gens(ph1, "*** after ph1.concatenate_assign(ph2) ***");

  return ok;
}

bool
test08() {
  Var A(0);
  Var B(1);

  Poly ph1(1);
  ph1.add_con(A >= 0);
  ph1.add_con(-A >= -2);
  ph1.minimize();
  ph1.add_gen(point(0*A));
  ph1.add_gen(point(2*A));

  Poly ph2(1, Spec_Elem::EMPTY);
  ph2.add_gen(point(10*A));

  print_cons(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");

  ph1.concatenate_assign(ph2);

  Poly known_result(2);
  known_result.add_con(A >= 0);
  known_result.add_con(-A >= -2);
  known_result.add_con(B == 10);

  bool ok = (ph1 == known_result);

  print_cons(ph1, "*** after ph1.concatenate_assign(ph2) ***");

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
END_MAIN
