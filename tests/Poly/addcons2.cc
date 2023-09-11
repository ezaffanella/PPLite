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
  Cons cs1(1, Con::zero_dim_false());
  Poly ph1;
  ph1.add_cons(cs1.begin(), cs1.end());

  Cons cs2;
  cs2.push_back(Linear_Expr() >= 7);
  Poly ph2;
  ph2.add_cons(cs2.begin(), cs2.end());

  Cons cs3;
  cs3.push_back(Linear_Expr() >= -3);
  Poly ph3;
  ph3.add_cons(cs3.begin(), cs3.end());

  Poly empty_result(0, Spec_Elem::EMPTY);
  Poly univ_result;

  bool ok = (ph1 == empty_result
             && ph2 == empty_result
             && ph3 == univ_result);

  return ok;
}

bool
test02() {
  Var A(0);
  Var B(1);

  Poly ph(2);
  ph.add_con(A >= 0);
  ph.add_con(B >= 0);

  Cons cs;
  Poly known_result(ph);
  ph.add_cons(cs);

  bool ok = (ph == known_result);
  return ok;
}

bool
test03() {
  Var x(0);
  Var y(1);

  Cons cs1;
  cs1.push_back(x + y >= 0);
  Poly ph(2);
  ph.add_cons(cs1);

  Cons cs2;
  cs2.push_back(Linear_Expr() == 1);
  ph.add_cons(cs2);

  Poly known_result(2, Spec_Elem::EMPTY);

  bool ok = (ph == known_result);
  return ok;
}

bool
test04() {
  Var x(0);

  Poly ph(3);
  ph.add_con(x >= 1);

  Cons cs;
  Poly computed_result(ph);
  computed_result.add_cons(cs);

  bool ok = (computed_result == ph);
  return ok;
}

bool
test05() {
  Var x(0);

  Poly ph(3, Spec_Elem::EMPTY);

  Cons cs;
  cs.push_back(x >= 4);

  ph.add_cons(cs);

  Poly computed_result(3, Spec_Elem::EMPTY);

  bool ok = (ph == computed_result);
  return ok;
}

bool
test06() {
  Var x(0);
  Var y(1);

  Poly ph(3);

  Cons cs;
  cs.push_back(x >= 4);
  cs.push_back(x - y >= 0);

  ph.add_cons(cs);

  Poly known_result(3);
  known_result.add_con(x >= 4);
  known_result.add_con(x - y >= 0);

  bool ok = (ph == known_result);
  return ok;
}

bool
test07() {
  Var x(0);
  Var y(1);

  Poly ph(3);
  ph.add_con(y >= 1);

  Cons cs;
  cs.push_back(x >= 0);
  cs.push_back(y <= 0);
  ph.add_cons(cs);

  Poly known_result(3, Spec_Elem::EMPTY);

  bool ok = (ph == known_result);
  return ok;
}

bool
test08() {
  Poly ph;
  ph.add_con(-2 >= Linear_Expr());

  Cons cs;
  cs.push_back(-1 >= Linear_Expr());

  ph.add_cons(cs);

  Poly known_result(0, Spec_Elem::EMPTY);

  bool ok = (known_result == ph);
  return ok;
}

bool
test09() {
  Var x(0);
  Var y(1);

  Gens gs;
  gs.push_back(point());
  gs.push_back(ray(x));
  gs.push_back(ray(x + y));

  Poly ph(2, Spec_Elem::EMPTY);
  ph.add_gens(gs);

  Cons cs;
  cs.push_back(x <= 3);
  ph.add_cons(cs);

  Poly known_result(2);
  known_result.add_con(y >= 0);
  known_result.add_con(x - y >= 0);
  known_result.add_con(x <= 3);

  bool ok = (known_result == ph);
  return ok;
}

bool
test10() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::EMPTY);

  Cons cs;
  cs.push_back(A - B >= 0);

  ph1.add_cons(cs);

  Poly known_result(2, Spec_Elem::EMPTY);

  bool ok = (ph1 == known_result);
  return ok;
}

bool
test11() {
  Var A(0);
  Var B(1);

  Poly ph(2);
  ph.copy_gens();
  ph.add_con(A >= 0);
  Poly copy_ph(ph);

  Cons cs1;
  cs1.push_back(A == 0);
  cs1.push_back(B >= 0);
  Cons cs2(cs1);

  ph.add_cons(cs1);
  copy_ph.add_cons(cs2);

  bool ok = (ph == copy_ph);
  return ok;
}

bool
test12() {
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

  Cons cs1 = ph2.copy_cons();
  Cons cs2 = ph2.copy_cons();

  ph1.add_cons(cs1);
  copy_ph1.add_cons(cs2);

  bool ok = (ph1 == copy_ph1);
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
