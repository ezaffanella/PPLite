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

#include <list>

bool
test01() {
  Var A(0);
  Poly ph1(1);
  ph1.add_con(A >= 0);
  Poly ph2(1);
  ph2.add_con(A >= -10);
  ph2.add_con(A >= -5);
  ph2.add_con(A >= -2);
  ph2.add_con(A >= -1);
  ph2.add_con(A >= 0);
  return ph1.is_necessarily_closed()
    && ph1.space_dim() == 1
    && ph1.equals(ph2);
}

bool
test02() {
  Var A(0);
  Var B(1);

  std::list<Con> cons;
  cons.push_back(A >= 0);
  cons.push_back(B == 5);

  Poly ph(2);
  ph.add_cons(cons.begin(), cons.end());

  Poly known_result(2);
  known_result.add_con(B == 5);
  known_result.add_con(A >= 0);

  bool ok = (ph == known_result);
  return ok;
}

bool
test03() {
  Var A(0);
  Var B(1);

  Poly ph(2, Spec_Elem::EMPTY);
  Poly known_result(ph);
  ph.add_con(A == B);
  bool ok = (ph == known_result);
  return ok;
}

bool
test04() {
  Var A(0);
  Var B(1);

  Poly ph(2);
  ph.add_con(A >= 0);
  ph.add_con(A <= 2);
  ph.add_con(A >= -1);
  ph.add_con(B >= 1);

  Poly known_result(2);
  known_result.add_con(A >= 0);
  known_result.add_con(A <= 2);
  known_result.add_con(B >= 1);

  bool ok = (ph == known_result);
  return ok;
}

bool
test05() {
  Var A(0);
  Var B(1);

  Poly ph(2);
  ph.add_con(B >= 0);
  ph.add_con(B <= 2);
  ph.add_con(A + B >= 1);
  ph.add_con(A - B >= -1);

  Poly known_result(2, Spec_Elem::EMPTY);
  known_result.add_gen(point(B));
  known_result.add_gen(ray(A));
  known_result.add_gen(point(A + 2*B));
  known_result.add_gen(point(A));

  bool ok = (ph == known_result);
  return ok;
}

bool
test06() {
  Var A(0);
  Var B(1);

  Poly ph(2);
  ph.add_con(B >= 0);
  ph.add_con(A >= 0);
  ph.add_con(A <= 2);

  Poly ph2(2, Spec_Elem::EMPTY);
  ph2.add_gen(point());
  ph2.add_gen(point(2*A));
  ph2.add_gen(line(B));

  bool ok = !(ph == ph2);
  return ok;
}

bool
test07() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::EMPTY);
  ph1.add_gen(point());
  ph1.add_gen(point(2*B));
  ph1.add_gen(ray(A));
  ph1.add_gen(ray(-A));
  ph1.minimize();

  Poly known_result(2, Spec_Elem::UNIVERSE);
  known_result.add_con(B <= 2);
  known_result.add_con(B >= 0);

  return ph1.is_necessarily_closed()
    && ph1 == known_result;
}

// Non maximal repr for efc would be a redundant ns
// which makes quick_equals check fail.
bool
test08() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::EMPTY, Topol::NNC);
  ph1.add_gen(point(A));
  ph1.add_gen(point(B));
  ph1.add_gen(closure_point());
  ph1.add_gen(ray(A));
  ph1.add_gen(ray(B));

  Poly ph2(2, Spec_Elem::UNIVERSE, Topol::NNC);
  ph2.add_con(A >= 0);
  ph2.add_con(A <= 1);
  ph2.minimize();

  ph1.intersection_assign(ph2);

  Poly known_result(2, Spec_Elem::UNIVERSE, Topol::NNC);
  known_result.add_con(A >= 0);
  known_result.add_con(B >= 0);
  known_result.add_con(A + B > 0);
  known_result.add_con(A <= 1);

  return ph1 == known_result;
}

bool
test09() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::UNIVERSE);
  ph1.add_con(B >= 0);
  ph1.add_con(A + B <= 2);
  ph1.add_con(A - B >= -2);
  ph1.minimize();
  ph1.add_con(A - B >= 2);

  Poly known_result(2, Spec_Elem::EMPTY);
  known_result.add_gen(point(2*A));

  return ph1 == known_result;
}

bool
test10() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::UNIVERSE);
  ph1.add_con(A >= 0);
  ph1.add_con(B >= 0);
  ph1.add_con(A <= 2);

  ph1.add_gen(ray(A));

  Poly known_result(2, Spec_Elem::UNIVERSE);
  known_result.add_con(A >= 0);
  known_result.add_con(B >= 0);

  return ph1 == known_result;
}

bool
test11() {
  Var A(0);

  Poly ph1(1, Topol::NNC);
  ph1.add_con(A > 0);
  ph1.add_con(A < 2);

  Poly ph2(1, Spec_Elem::EMPTY, Topol::NNC);
  ph2.add_gens(ph1.copy_gens());
  ph2.ascii_dump(nout);

  bool ok = (!ph2.is_empty());
  return ok;
}

bool
test12() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Topol::NNC);
  ph1.add_con(A > 0);
  ph1.add_con(A < 1);
  ph1.add_gen(ray(A));

  ph1.add_con(A > 2);
  ph1.minimize();
  ph1.ascii_dump(nout);

  bool ok = (ph1.check_inv());
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

