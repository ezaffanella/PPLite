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

  Poly ph1(2, Topol::NNC);
  ph1.add_con(A >= 0);
  ph1.add_con(B >= 0);
  ph1.add_con(A + B > 0);

  ph1.add_gen(point());

  Poly known_result(2, Topol::NNC);
  known_result.add_con(A >= 0);
  known_result.add_con(B >= 0);

  ph1.minimize();
  known_result.minimize();

  return ph1 == known_result;
}


bool
test02() {
  Var A(0);
  Var B(1);
  Var C(2);

  Poly ph1(2, Topol::NNC);
  ph1.add_con(A >= 0);
  ph1.add_con(B >= 0);
  ph1.add_con(B <= 2);
  ph1.add_con(A + B > 0);

  ph1.add_gen(point());

  Poly known_result(2, Topol::NNC);
  known_result.add_con(A >= 0);
  known_result.add_con(B >= 0);
  known_result.add_con(B <= 2);

  return ph1 == known_result;
}

bool
test03() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Topol::NNC);
  ph1.add_con(A >= 0);
  ph1.add_con(B >= 0);
  ph1.add_con(A + B > 0);

  ph1.add_gen(ray(-A));

  Poly known_result(2, Topol::NNC);
  known_result.add_con(B >= 0);

  return ph1 == known_result;
}

bool
test04() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Topol::NNC);
  ph1.add_con(A >= 0);
  ph1.add_con(B >= 0);
  ph1.add_con(A + B > 0);

  ph1.add_gen(line(-A));

  Poly known_result(2, Topol::NNC);
  known_result.add_con(B >= 0);

  return ph1 == known_result;
}

bool
test05() {
  Var A(0);
  Var B(1);
  Var C(2);

  Poly ph1(3, Topol::NNC);
  ph1.add_con(A >= 0);
  ph1.add_con(B >= 0);
  ph1.add_con(A + B > 0);

  ph1.add_gen(ray(-A - B));

  Poly known_result(3, Topol::NNC);

  return ph1 == known_result;
}

bool
test06() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Topol::NNC);
  ph1.add_con(A >= 0);
  ph1.add_con(B >= 0);

  ph1.add_gen(ray(A));
  ph1.add_gen(closure_point());

  Poly known_result(2, Topol::NNC);
  known_result.add_con(A >= 0);
  known_result.add_con(B >= 0);

  return ph1 == known_result;
}

bool
test07() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Topol::NNC);
  ph1.add_con(A >= 1);
  ph1.add_con(B >= 0);

  ph1.add_gen(closure_point());

  Poly known_result(2, Topol::NNC);
  known_result.add_con(A > 0);
  known_result.add_con(B >= 0);

  return ph1 == known_result;
}

bool
test08() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Topol::NNC);
  ph1.add_con(A >= 1);
  ph1.add_con(B >= 0);

  ph1.add_gen(closure_point());
  ph1.add_gen(point(B));

  Poly known_result(2, Topol::NNC);
  known_result.add_con(A >= 0);
  known_result.add_con(B >= 0);
  known_result.add_con(A + B > 0);

  return ph1 == known_result;
}

bool
test09() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Topol::NNC);
  ph1.add_con(A > 0);
  ph1.add_con(B > 0);
  ph1.add_con(A + B < 1);

  ph1.add_gen(point());
  ph1.add_gen(point(A));
  ph1.add_gen(point(B));

  Poly known_result(2, Topol::NNC);
  known_result.add_con(A >= 0);
  known_result.add_con(B >= 0);
  known_result.add_con(A + B <= 1);

  return ph1 == known_result;
}

bool
test10() {
  Var A(0);
  Var B(1);
  Var C(2);

  Poly ph1(3, Topol::NNC);
  ph1.add_con(A >= 0);
  ph1.add_con(B >= 0);
  ph1.add_con(A + B > 0);

  ph1.add_gen(ray(-A - B));

  Poly known_result(3, Spec_Elem::EMPTY, Topol::NNC);
  known_result.add_gen(point());
  known_result.add_gen(line(A));
  known_result.add_gen(line(B));
  known_result.add_gen(line(C));

  return ph1 == known_result;
}

bool
test11() {
  Var A(0);
  Var B(1);
  Var C(2);

  Poly ph1(3, Topol::NNC);
  ph1.add_con(C == 0);
  ph1.add_con(A >= 0);
  ph1.add_con(B >= 0);
  ph1.add_con(A + B > 0);

  ph1.add_gen(ray(-A - B));

  Poly known_result(3, Spec_Elem::EMPTY, Topol::NNC);
  known_result.add_gen(point());
  known_result.add_gen(line(A));
  known_result.add_gen(line(B));

  return ph1 == known_result;
}

bool
test12() {
  Var A(0);
  Var B(1);

  Poly ph1(5, Topol::NNC);
  ph1.add_con(A >= 0);
  ph1.add_con(B >= 0);
  ph1.add_con(A <= 2);
  ph1.add_con(B <= 3);

  ph1.add_gen(line(A));

  Poly known_result(5, Topol::NNC);
  known_result.add_con(B <= 3);
  known_result.add_con(B >= 0);

  return ph1 == known_result;
}

bool
test13() {
  Var A(0);
  Var B(1);

  Poly ph1(5, Topol::NNC);
  ph1.add_con(A >= 0);
  ph1.add_con(B >= 0);
  ph1.add_con(A + B < 3);

  ph1.add_gen(ray(-A -B));

  Poly known_result(5, Spec_Elem::EMPTY, Topol::NNC);
  known_result.add_gen(point());
  known_result.add_gen(closure_point(3*A));
  known_result.add_gen(closure_point(3*B));
  known_result.add_gen(ray(-A -B));
  known_result.add_gen(line(Var(2)));
  known_result.add_gen(line(Var(3)));
  known_result.add_gen(line(Var(4)));

  return ph1 == known_result;
}


bool
test14() {
  Var A(0);
  Var B(1);

  Poly ph1(5, Topol::NNC);
  ph1.add_con(A >= 0);
  ph1.add_con(B >= 0);
  ph1.add_con(A + B <= 4);
  ph1.add_con(A < 4);
  ph1.add_con(B < 4);

  ph1.add_gen(ray(-A -B));

  Poly known_result(5, Spec_Elem::EMPTY, Topol::NNC);
  known_result.add_gen(point(2*A + 2*B));
  known_result.add_gen(closure_point(4*A));
  known_result.add_gen(closure_point(4*B));
  known_result.add_gen(ray(-A -B));
  known_result.add_gen(line(Var(2)));
  known_result.add_gen(line(Var(3)));
  known_result.add_gen(line(Var(4)));

  return ph1 == known_result;
}


bool
test15() {
  Var A(0);
  Var B(1);

  Poly ph1(5, Topol::NNC);
  ph1.add_con(A >= 0);
  ph1.add_con(B >= 0);
  ph1.add_con(A <= 3);
  ph1.add_con(B <= 3);
  ph1.add_con(A - B > -2);
  ph1.add_con(A - B <  2);
  ph1.add_con(A + B >= 1);
  ph1.add_con(A + B <= 5);

  ph1.add_gen(point(2*A));

  Poly known_result(5, Spec_Elem::EMPTY, Topol::NNC);
  known_result.add_gen(point(A));
  known_result.add_gen(point(2*A));
  known_result.add_gen(closure_point(3*A + B));
  known_result.add_gen(point(3*A + 2*B));
  known_result.add_gen(point(2*A + 3*B));
  known_result.add_gen(closure_point(A + 3*B));
  known_result.add_gen(closure_point(2*B));
  known_result.add_gen(point(B));
  known_result.add_gen(line(Var(2)));
  known_result.add_gen(line(Var(3)));
  known_result.add_gen(line(Var(4)));

  return ph1 == known_result;
}

bool
test16() {
  Var A(0);
  Var B(1);

  Poly ph1(5, Spec_Elem::EMPTY, Topol::NNC);
  ph1.add_gen(point(A));
  ph1.add_gen(closure_point(2*A));
  ph1.add_gen(closure_point(3*A + B));
  ph1.add_gen(point(3*A + 2*B));
  ph1.add_gen(point(2*A + 3*B));
  ph1.add_gen(point(A + 3*B));
  ph1.add_gen(point(2*B));
  ph1.add_gen(point(B));
  ph1.add_gen(line(Var(2)));
  ph1.add_gen(line(Var(3)));
  ph1.add_gen(line(Var(4)));

  ph1.add_con(A - B > 1);

  Poly known_result(5, Topol::NNC);
  known_result.add_con(A <= 3);
  known_result.add_con(B >= 0);
  known_result.add_con(A - B < 2);
  known_result.add_con(A - B > 1);

  return ph1 == known_result;
}

bool
test17() {
  Var A(0);
  Var B(1);

  Poly ph1(5, Spec_Elem::EMPTY, Topol::NNC);
  ph1.add_gen(point(A));
  ph1.add_gen(closure_point(2*A));
  ph1.add_gen(closure_point(3*A + B));
  ph1.add_gen(point(3*A + 2*B));
  ph1.add_gen(point(2*A + 3*B));
  ph1.add_gen(point(A + 3*B));
  ph1.add_gen(point(2*B));
  ph1.add_gen(point(B));
  ph1.add_gen(line(Var(2)));
  ph1.add_gen(line(Var(3)));
  ph1.add_gen(line(Var(4)));

  ph1.add_con(A - B > 0);

  Poly known_result(5, Topol::NNC);
  known_result.add_con(A <= 3);
  known_result.add_con(B >= 0);
  known_result.add_con(A + B >= 1);
  known_result.add_con(A + B <= 5);
  known_result.add_con(A - B < 2);
  known_result.add_con(A - B > 0);

  return ph1 == known_result;
}

bool
test18() {
  Var A(0);
  Var B(1);

  Poly ph1(5, Topol::NNC);
  ph1.add_con(A >= 0);
  ph1.add_con(A <= 5);
  ph1.add_con(B >= 0);
  ph1.add_con(B <= 5);

  ph1.add_gen(ray(A));

  Poly known_result(5, Topol::NNC);
  known_result.add_con(A >= 0);
  known_result.add_con(B >= 0);
  known_result.add_con(B <= 5);

  return ph1 == known_result;
}

bool
test19() {
  Var A(0);
  Var B(1);

  Poly ph1(5, Topol::NNC);
  ph1.add_con(A > 0);
  ph1.add_con(A < 5);
  ph1.add_con(B >= 0);
  ph1.add_con(B <= 4);

  ph1.add_gen(ray(A));
  ph1.add_gen(closure_point(5*B));

  ph1.minimize();

  Poly known_result(5, Topol::NNC);
  known_result.add_con(A > 0);
  known_result.add_con(B >= 0);
  known_result.add_con(B < 5);

  return ph1 == known_result;
}

bool
test20() {
  Var A(0);
  Var B(1);

  Poly ph1(5, Topol::NNC);
  ph1.add_con(A > 0);
  ph1.add_con(A < 5);
  ph1.add_con(B >= 0);
  ph1.add_con(B <= 4);

  ph1.add_gen(line(A));
  ph1.add_gen(closure_point(5*B));

  Poly known_result(5, Topol::NNC);
  known_result.add_con(B >= 0);
  known_result.add_con(B < 5);

  return ph1 == known_result;
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
