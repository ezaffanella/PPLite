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
  ph1.add_con(x >= 1);
  ph1.add_con(x <= 3);
  ph1.add_con(y >= 1);
  ph1.add_con(y <= 3);

  Poly ph2(2);
  ph2.add_con(y == 5);

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  ph1.time_elapse_assign(ph2);

  Poly known_result(2, Spec_Elem::EMPTY);
  known_result.add_gen(point(y));
  known_result.add_gen(ray(y));
  known_result.add_gen(line(x));

  bool ok = (ph1 == known_result);

  print_gens(ph1, "*** ph1.time_elapse_assign(ph2) ***");

  return ok;
}

bool
test02() {
  Poly ph1(0, Spec_Elem::EMPTY);
  Poly ph2;

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  ph1.time_elapse_assign(ph2);

  Poly ph3(2, Spec_Elem::EMPTY);
  Poly ph4(2);

  print_cons(ph3, "*** ph3 ***");
  print_cons(ph4, "*** ph4 ***");

  ph3.time_elapse_assign(ph4);

  Poly ph5(2);
  Poly ph6(2, Spec_Elem::EMPTY);

  print_cons(ph5, "*** ph5 ***");
  print_cons(ph6, "*** ph6 ***");

  ph5.time_elapse_assign(ph6);

  bool ok = ph1.is_empty() && ph3.is_empty() && ph5.is_empty();

  print_gens(ph1, "*** ph1.time_elapse_assign(ph2) ***");
  print_gens(ph3, "*** ph3.time_elapse_assign(ph4) ***");
  print_gens(ph5, "*** ph5.time_elapse_assign(ph6) ***");

  return ok;
}

bool
test03() {
  Var A(0);
  Var B(1);

  Gens gs1;
  gs1.push_back(point(A));
  gs1.push_back(point(A + B));
  Poly ph1(2, Spec_Elem::EMPTY);
  ph1.add_gens(gs1);

  Gens gs2;
  gs2.push_back(point(0*B));
  Poly ph2(2, Spec_Elem::EMPTY);
  ph2.add_gens(gs2);

  print_gens(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");

  Poly known_result(ph1);

  ph1.time_elapse_assign(ph2);

  bool ok = (ph1 == known_result);

  print_gens(ph1, "*** ph1 ***");

  return ok;
}

bool
test04() {
  Poly ph1;
  Poly ph2(0, Spec_Elem::EMPTY);

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  ph1.time_elapse_assign(ph2);

  Poly known_result(0, Spec_Elem::EMPTY);

  bool ok = (ph1 == known_result);

  print_cons(ph1, "*** after ph1.time_elapse_assign(ph2) ***");

  return ok;
}

bool
test05() {
  Var A(0);
  Var B(1);

  Poly ph1(2);
  ph1.minimize();
  ph1.add_con(A == 0);
  ph1.add_con(B == 0);

  Poly ph2(2);
  ph2.minimize();
  ph2.add_con(A == 2);
  ph2.add_con(B == 2);

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  ph1.time_elapse_assign(ph2);

  Poly ph3(2, Spec_Elem::EMPTY);
  ph3.add_gen(point());
  ph3.minimize();

  Poly ph4(2, Spec_Elem::EMPTY);
  ph4.add_gen(point(2*A + 2*B));

  print_gens(ph3, "*** ph3 ***");
  print_gens(ph4, "*** ph4 ***");

  ph3.time_elapse_assign(ph4);

  bool ok = (ph1 ==  ph3);

  print_gens(ph1, "*** after ph1.time_elapse_assign(ph2) ***");
  print_gens(ph3, "*** after ph3.time_elapse_assign(ph4) ***");

  return ok;
}

bool
test06() {
  Var x(0);
  Var y(1);

  Poly ph1(Topol::NNC, 2);
  ph1.add_con(x >= 0);
  ph1.add_con(y >= 0);
  ph1.add_con(x + y <= 2);

  Poly ph2(Topol::NNC, 2);
  ph2.add_con(x > 2);
  ph2.add_con(x < 4);
  ph2.add_con(y == 3);

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  ph1.time_elapse_assign(ph2);

  Gens known_gs;
  known_gs.push_back(point());
  known_gs.push_back(point(2*x));
  known_gs.push_back(point(2*y));
  known_gs.push_back(ray(2*x + 3*y));
  known_gs.push_back(ray(4*x + 3*y));

  Poly known_result(Topol::NNC, Spec_Elem::EMPTY, 2);
  known_result.add_gens(known_gs);

  bool ok = (ph1 == known_result);

  print_gens(ph1, "*** ph1.time_elapse_assign(ph2) ***");

  return ok;
}

bool
test07() {
  Var x(0);
  Var y(1);

  Cons cs1;
  cs1.push_back(x > 3);
  cs1.push_back(y > 3);
  Poly ph(Topol::NNC, 2);
  ph.add_cons(cs1);

  Poly ph1(ph);

  Gens gs;
  gs.push_back(point(x + y));
  Poly ph2(Topol::NNC, Spec_Elem::EMPTY, 2);
  ph2.add_gens(gs);

  print_cons(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");

  ph1.time_elapse_assign(ph2);

  bool ok = (ph1 == ph);

  print_gens(ph1, "*** ph1.time_elapse_assign(ph2) ***");

 return ok;
}

bool
test08() {
  Var x(0);
  Var y(1);

  Poly ph1(Topol::NNC, 2);
  ph1.add_con(x == 3);
  ph1.add_con(y > 2);

  Poly ph2(Topol::NNC, 2);
  ph2.add_con(x > 3);
  ph2.add_con(y > 2);

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  ph1.time_elapse_assign(ph2);

  Gens gs;
  gs.push_back(closure_point(3*x + 2*y));
  gs.push_back(point(3*x + 3*y));
  gs.push_back(ray(y));
  gs.push_back(ray(x));

  Poly known_result(Topol::NNC, Spec_Elem::EMPTY, 2);
  known_result.add_gens(gs);

  bool ok = (ph1 == known_result);

  print_gens(ph1, "*** ph1.time_elapse_assign(ph2) ***");

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
