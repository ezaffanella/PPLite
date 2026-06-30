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
  Var x(0);

  Poly ph1(1, Spec_Elem::UNIVERSE, Topol::NNC);
  ph1.add_con(x > 0);
  ph1.add_con(x <= 1);

  print_cons(ph1, "*** ph1 ***");

  Poly ph2(1, Spec_Elem::UNIVERSE, Topol::NNC);
  ph2.add_con(x == 0);

  print_cons(ph2, "*** ph2 ***");

  Poly kr(1, Spec_Elem::UNIVERSE, Topol::NNC);
  kr.add_con(x >= 0);
  kr.add_con(x <= 1);

  bool ok = ph1.join_assign_if_exact(ph2);
  ok &= (ph1 == kr);

  print_cons(ph1, "*** ph1.join_assign_if_exact(ph2) ***");

  return ok;
}

bool
test02() {
  Var x(0);
  Var y(1);

  Poly ph1(2, Spec_Elem::UNIVERSE, Topol::NNC);
  ph1.add_con(x > 0);
  ph1.add_con(x <= 2);
  ph1.add_con(y >= 0);
  ph1.add_con(y <= 2);

  print_cons(ph1, "*** ph1 ***");

  Poly ph2(2, Spec_Elem::UNIVERSE, Topol::NNC);
  ph2.add_con(x >= 0);
  ph2.add_con(x <= 2);
  ph2.add_con(y == 1);

  print_cons(ph2, "*** ph2 ***");

  Poly kr(ph1);

  bool ok = !ph1.join_assign_if_exact(ph2);
  ok &= (ph1 == kr);

  print_cons(ph1, "*** ph1.join_assign_if_exact(ph2) ***");

  return ok;
}

bool
test03() {
  Var x(0);
  Var y(1);

  Poly ph1(2, Spec_Elem::UNIVERSE, Topol::NNC);
  ph1.add_con(x >= 0);
  ph1.add_con(x <= 2);
  ph1.add_con(y >= 0);
  ph1.add_con(y <= 2);
  ph1.add_con(x + y > 0);
  ph1.add_con(x + y < 4);

  print_cons(ph1, "*** ph1 ***");

  Poly ph2(2, Spec_Elem::UNIVERSE, Topol::NNC);
  ph2.add_con(x == y);
  ph2.add_con(x <= 2);
  ph2.add_con(y >= 0);

  print_cons(ph2, "*** ph2 ***");

  Poly kr(2, Spec_Elem::UNIVERSE, Topol::NNC);
  kr.add_con(x >= 0);
  kr.add_con(x <= 2);
  kr.add_con(y >= 0);
  kr.add_con(y <= 2);

  bool ok = ph1.join_assign_if_exact(ph2);
  ok &= (ph1 == kr);

  print_cons(ph1, "*** ph1.join_assign_if_exact(ph2) ***");

  return ok;
}

bool
test04() {
  Var x(0);
  Var y(1);

  Poly ph1(2, Spec_Elem::UNIVERSE, Topol::NNC);
  ph1.add_con(x - y >= 0);
  ph1.add_con(x + y >= 0);
  ph1.add_con(4*x < 1);

  print_cons(ph1, "*** ph1 ***");

  Poly ph2(2, Spec_Elem::UNIVERSE, Topol::NNC);
  ph2.add_con(x - y >= 0);
  ph2.add_con(x + y >= 0);
  ph2.add_con(4*x > 1);

  print_cons(ph2, "*** ph2 ***");

  Poly kr(ph1);

  bool ok = !ph1.join_assign_if_exact(ph2);
  ok &= (ph1 == kr);

  print_cons(ph1, "*** ph1.join_assign_if_exact(ph2) ***");

  return ok;
}

bool
test05() {
  Var x(0);
  Var y(1);

  Poly ph1(2, Spec_Elem::UNIVERSE, Topol::NNC);
  ph1.add_con(x > 0);
  ph1.add_con(y > 0);
  ph1.add_con(x < 2);
  ph1.add_con(y < 2);

  print_cons(ph1, "*** ph1 ***");

  Poly ph2(2, Spec_Elem::UNIVERSE, Topol::NNC);
  ph2.add_con(x == 2);
  ph2.add_con(y == 1);

  print_cons(ph2, "*** ph2 ***");

  Poly kr(ph1);

  bool ok = !ph1.join_assign_if_exact(ph2);
  ok &= (ph1 == kr);

  print_cons(ph1, "*** ph1.join_assign_if_exact(ph2) ***");

  return ok;
}

bool
test06() {
  Var x(0);
  Var y(1);

  Poly ph1(2, Spec_Elem::UNIVERSE, Topol::NNC);
  ph1.add_con(x > 0);
  ph1.add_con(y > 0);
  ph1.add_con(x < 2);
  ph1.add_con(y < 2);

  print_cons(ph1, "*** ph1 ***");

  Poly ph2(2, Spec_Elem::UNIVERSE, Topol::NNC);
  ph2.add_con(x == 2);
  ph2.add_con(y > 0);
  ph2.add_con(y < 2);

  print_cons(ph2, "*** ph2 ***");

  Poly kr(2, Spec_Elem::UNIVERSE, Topol::NNC);
  kr.add_con(x > 0);
  kr.add_con(y > 0);
  kr.add_con(x <= 2);
  kr.add_con(y < 2);

  bool ok = ph1.join_assign_if_exact(ph2);
  ok &= (ph1 == kr);

  print_cons(ph1, "*** ph1.join_assign_if_exact(ph2) ***");

  return ok;
}

bool
test07() {
  Var x(0);
  Var y(1);

  Poly ph1(2, Spec_Elem::UNIVERSE, Topol::NNC);
  ph1.add_con(x >= 0);
  ph1.add_con(x < 2);

  print_cons(ph1, "*** ph1 ***");

  Poly ph2(2, Spec_Elem::UNIVERSE, Topol::NNC);
  ph2.add_con(x == 2);
  ph2.add_con(y == 0);

  print_cons(ph2, "*** ph2 ***");

  Poly kr(ph1);

  bool ok = !ph1.join_assign_if_exact(ph2);
  ok &= (ph1 == kr);

  print_cons(ph1, "*** ph1.join_assign_if_exact(ph2) ***");

  return ok;
}

bool
test08() {
  Var x(0);
  Var y(1);

  Poly ph1(2, Spec_Elem::UNIVERSE, Topol::NNC);
  ph1.add_con(x >= 0);
  ph1.add_con(x < 1);
  ph1.add_con(y >= 0);
  ph1.add_con(y <= 2);

  print_cons(ph1, "*** ph1 ***");

  Poly ph2(2, Spec_Elem::UNIVERSE, Topol::NNC);
  ph2.add_con(x > 1);
  ph2.add_con(x <= 2);
  ph2.add_con(y >= 0);
  ph2.add_con(y <= 2);

  print_cons(ph2, "*** ph2 ***");

  Poly kr(ph1);

  bool ok = !ph1.join_assign_if_exact(ph2);
  ok &= (ph1 == kr);

  print_cons(ph1, "*** ph1.join_assign_if_exact(ph2) ***");

  return ok;
}

bool
test09() {
  Var x(0);
  Var y(1);

  Poly ph1(2, Spec_Elem::EMPTY, Topol::NNC);
  ph1.add_gen(point(0*x + 0*y));
  ph1.add_gen(point(0*x + 1*y));
  ph1.add_gen(closure_point(2*x + 2*y));
  ph1.add_gen(closure_point(3*x + 0*y));

  print_cons(ph1, "*** ph1 ***");

  Poly ph2(2, Spec_Elem::EMPTY, Topol::NNC);
  ph2.add_gen(point(4*x + 0*y));
  ph2.add_gen(point(4*x + 1*y));
  ph2.add_gen(closure_point(2*x + 2*y));
  ph2.add_gen(closure_point(1*x + 0*y));

  print_cons(ph2, "*** ph2 ***");

  Poly kr(ph1);
  kr.add_gen(point(0*x + 0*y));
  kr.add_gen(point(0*x + 1*y));
  kr.add_gen(point(4*x + 0*y));
  kr.add_gen(point(4*x + 1*y));
  kr.add_gen(closure_point(2*x + 2*y));

  bool ok = ph1.join_assign_if_exact(ph2);
  ok &= (ph1 == kr);

  print_cons(ph1, "*** ph1.join_assign_if_exact(ph2) ***");

  return ok;
}

bool
test10() {
  Var x(0);
  Var y(1);

  Poly ph1(2, Spec_Elem::UNIVERSE, Topol::NNC);
  ph1.add_con(x > 0);
  ph1.add_con(x < 1);
  ph1.add_con(y > 0);
  ph1.add_con(y < 2);

  print_cons(ph1, "*** ph1 ***");

  Poly ph2(2, Spec_Elem::UNIVERSE, Topol::NNC);
  ph2.add_con(x > 1);
  ph2.add_con(x < 2);
  ph2.add_con(y > 0);
  ph2.add_con(y < 2);

  print_cons(ph2, "*** ph2 ***");

  Poly kr(ph1);

  bool ok = !ph1.join_assign_if_exact(ph2);
  ok &= (ph1 == kr);

  print_cons(ph1, "*** ph1.join_assign_if_exact(ph2) ***");

  return ok;
}

bool
test11() {
  Var x(0);
  Var y(1);

  Poly ph1(2, Spec_Elem::UNIVERSE, Topol::NNC);
  ph1.add_con(x > 0);
  ph1.add_con(x < 1);
  ph1.add_con(y > 0);
  ph1.add_con(y < 2);

  print_cons(ph1, "*** ph1 ***");

  Poly ph2(2, Spec_Elem::UNIVERSE, Topol::NNC);
  ph2.add_con(x >= 1);
  ph2.add_con(x < 2);
  ph2.add_con(y > 0);
  ph2.add_con(y < 2);

  print_cons(ph2, "*** ph2 ***");

  Poly kr(2, Spec_Elem::UNIVERSE, Topol::NNC);
  kr.add_con(x > 0);
  kr.add_con(x < 2);
  kr.add_con(y > 0);
  kr.add_con(y < 2);

  bool ok = ph1.join_assign_if_exact(ph2);
  ok &= (ph1 == kr);

  print_cons(ph1, "*** ph1.join_assign_if_exact(ph2) ***");

  return ok;
}

bool
test12() {
  Var x(0);
  Var y(1);
  Var z(2);

  Poly ph1(3, Spec_Elem::UNIVERSE, Topol::NNC);
  ph1.add_con(x > 0);
  ph1.add_con(x <= 1);
  ph1.add_con(y > 0);
  ph1.add_con(y < 2);
  ph1.add_con(z > 0);
  ph1.add_con(z < 2);

  print_cons(ph1, "*** ph1 ***");

  Poly ph2(3, Spec_Elem::UNIVERSE, Topol::NNC);
  ph2.add_con(x >= 1);
  ph2.add_con(x < 2);
  ph2.add_con(y > 0);
  ph2.add_con(y < 2);
  ph2.add_con(z > 0);
  ph2.add_con(z < 2);

  print_cons(ph2, "*** ph2 ***");

  Poly kr(3, Spec_Elem::UNIVERSE, Topol::NNC);
  kr.add_con(x > 0);
  kr.add_con(x < 2);
  kr.add_con(y > 0);
  kr.add_con(y < 2);
  kr.add_con(z > 0);
  kr.add_con(z < 2);

  bool ok = ph1.join_assign_if_exact(ph2);
  ok &= (ph1 == kr);

  print_cons(ph1, "*** ph1.join_assign_if_exact(ph2) ***");

  return ok;
}

bool
test13() {
  Var x(0);
  Var y(1);
  Var z(2);

  Poly ph1(3, Spec_Elem::UNIVERSE, Topol::NNC);
  ph1.add_con(x > 0);
  ph1.add_con(x <= 1);
  ph1.add_con(y > 0);
  ph1.add_con(y < 2);
  ph1.add_con(z > 0);
  ph1.add_con(z <= 2);
  ph1.add_con(x + z < 3);

  ph1.minimize();
  print_cons(ph1, "*** ph1 ***");
  print_gens(ph1, "*** ph1 ***");

  Poly ph2(3, Spec_Elem::UNIVERSE, Topol::NNC);
  ph2.add_con(x >= 1);
  ph2.add_con(x < 2);
  ph2.add_con(y > 0);
  ph2.add_con(y < 2);
  ph2.add_con(z > 0);
  ph2.add_con(z <= 2);
  ph1.add_con(x - z > -1);

  ph2.minimize();
  print_cons(ph2, "*** ph2 ***");
  print_gens(ph2, "*** ph2 ***");

  Poly kr(ph1);

  bool ok = !ph1.join_assign_if_exact(ph2);
  ok &= (ph1 == kr);

  print_cons(ph1, "*** ph1.join_assign_if_exact(ph2) ***");

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
END_MAIN
