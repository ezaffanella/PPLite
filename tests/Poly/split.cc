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
  Var B(1);
  Cons cs;
  cs.push_back(A >= 0);
  cs.push_back(B >= 0);
  cs.push_back(A <= 4);
  cs.push_back(B <= 4);

  Poly ph1(2, Topol::NNC);
  ph1.add_cons(cs);

  Con c(A > 4);

  Poly ph2 = ph1.split(c);

  Poly knr1(2, Topol::NNC);
  knr1.add_cons(cs);
  knr1.add_con(A > 4);
  knr1.minimize();

  Poly knr2(2, Topol::NNC);
  knr2.add_cons(cs);
  knr2.add_con(A <= 4);
  knr2.minimize();

  return (ph1 == knr1 && ph2 == knr2);
}

bool
test02() {
  Var A(0);
  Var B(1);
  Cons cs;
  cs.push_back(A >= 0);
  cs.push_back(B >= 0);
  cs.push_back(A <= 4);
  cs.push_back(B <= 4);

  Poly ph1(2, Topol::NNC);
  ph1.add_cons(cs);

  Con c(A > 2);

  Poly ph2 = ph1.split(c);

  Poly knr1(2, Topol::NNC);
  knr1.add_cons(cs);
  knr1.add_con(A > 2);
  knr1.minimize();

  Poly knr2(2, Topol::NNC);
  knr2.add_cons(cs);
  knr2.add_con(A <= 2);
  knr2.minimize();

  return (ph1 == knr1 && ph2 == knr2);
}

bool
test03() {
  Var A(0);
  Var B(1);
  Cons cs;
  cs.push_back(A >= 0);
  cs.push_back(B >= 0);
  cs.push_back(A <= 4);
  cs.push_back(B <= 4);

  Poly ph1(2, Topol::NNC);
  ph1.add_cons(cs);

  Con c(A <= 0);

  Poly ph2 = ph1.split(c);

  Poly knr1(2, Topol::NNC);
  knr1.add_cons(cs);
  knr1.add_con(A <= 0);
  knr1.minimize();

  Poly knr2(2, Topol::NNC);
  knr2.add_cons(cs);
  knr2.add_con(A > 0);
  knr2.minimize();

  return (ph1 == knr1 && ph2 == knr2);
}

bool
test04() {
  Var A(0);
  Var B(1);
  Cons cs;
  cs.push_back(A >= 0);
  cs.push_back(B >= 0);
  cs.push_back(A <= 4);
  cs.push_back(B <= 4);

  Poly ph1(2, Topol::NNC);
  ph1.add_cons(cs);

  Con c(A <= 2);

  Poly ph2 = ph1.split(c);

  Poly knr1(2, Topol::NNC);
  knr1.add_cons(cs);
  knr1.add_con(A <= 2);
  knr1.minimize();

  Poly knr2(2, Topol::NNC);
  knr2.add_cons(cs);
  knr2.add_con(A > 2);
  knr2.minimize();

  return (ph1 == knr1 && ph2 == knr2);
}

bool
test05() {
  Var A(0);
  Var B(1);
  Cons cs;
  cs.push_back(A >= 0);
  cs.push_back(B >= 0);
  cs.push_back(A <= 4);
  cs.push_back(B <= 4);

  Poly ph1(2);
  ph1.add_cons(cs);

  Con c(A >= 4);

  Poly ph2 = ph1.split(c);

  Poly knr1(2);
  knr1.add_cons(cs);
  knr1.add_con(A >= 4);
  knr1.minimize();

  Poly knr2(2);
  knr2.add_cons(cs);
  knr2.add_con(A <= 4);
  knr2.minimize();

  return (ph1 == knr1 && ph2 == knr2);
}

bool
test06() {
  Var A(0);
  Var B(1);
  Cons cs;
  cs.push_back(A >= 0);
  cs.push_back(B >= 0);
  cs.push_back(A <= 4);
  cs.push_back(B <= 4);

  Poly ph1(2);
  ph1.add_cons(cs);

  Con c(A >= 2);

  Poly ph2 = ph1.split(c);

  Poly knr1(2);
  knr1.add_cons(cs);
  knr1.add_con(A >= 2);
  knr1.minimize();

  Poly knr2(2);
  knr2.add_cons(cs);
  knr2.add_con(A <= 2);
  knr2.minimize();

  return (ph1 == knr1 && ph2 == knr2);
}

bool
test07() {
  Var A(0);
  Var B(1);
  Cons cs;
  cs.push_back(A >= 0);
  cs.push_back(B >= 0);
  cs.push_back(A <= 4);
  cs.push_back(B <= 4);

  Poly ph1(2);
  ph1.add_cons(cs);

  Con c(A <= 0);

  Poly ph2 = ph1.split(c);

  Poly knr1(2);
  knr1.add_cons(cs);
  knr1.add_con(A <= 0);
  knr1.minimize();

  Poly knr2(2);
  knr2.add_cons(cs);
  knr2.add_con(A >= 0);
  knr2.minimize();

  return (ph1 == knr1 && ph2 == knr2);
}

bool
test08() {
  Var A(0);
  Var B(1);
  Cons cs;
  cs.push_back(A >= 0);
  cs.push_back(B >= 0);
  cs.push_back(A <= 4);
  cs.push_back(B <= 4);

  Poly ph1(2);
  ph1.add_cons(cs);

  Con c(A <= 2);

  Poly ph2 = ph1.split(c);

  Poly knr1(2);
  knr1.add_cons(cs);
  knr1.add_con(A <= 2);
  knr1.minimize();

  Poly knr2(2);
  knr2.add_cons(cs);
  knr2.add_con(A >= 2);
  knr2.minimize();

  return (ph1 == knr1 && ph2 == knr2);
}

bool
test09() {
  Var A(0);
  Var B(1);
  Cons cs;
  cs.push_back(A >= 0);
  cs.push_back(B > 0);
  cs.push_back(A <= 4);
  cs.push_back(B < 4);

  Poly ph1(2, Topol::NNC);
  ph1.add_cons(cs);

  Con c(B <= 2);

  Poly ph2 = ph1.split(c);

  Poly knr1(2, Topol::NNC);
  knr1.add_cons(cs);
  knr1.add_con(B <= 2);
  knr1.minimize();

  Poly knr2(2, Topol::NNC);
  knr2.add_cons(cs);
  knr2.add_con(B > 2);
  knr2.minimize();

  return (ph1 == knr1 && ph2 == knr2);
}

bool
test10() {
  Var A(0);
  Var B(1);
  Cons cs;
  cs.push_back(A > 0);
  cs.push_back(B > 0);
  cs.push_back(A < 4);
  cs.push_back(B < 4);

  Poly ph1(2, Topol::NNC);
  ph1.add_cons(cs);

  Con c(B < 2);

  Poly ph2 = ph1.split(c);

  Poly knr1(2, Topol::NNC);
  knr1.add_cons(cs);
  knr1.add_con(B < 2);
  knr1.minimize();

  Poly knr2(2, Topol::NNC);
  knr2.add_cons(cs);
  knr2.add_con(B >= 2);
  knr2.minimize();

  return (ph1 == knr1 && ph2 == knr2);
}

bool
test11() {
  Var A(0);
  Var B(1);
  Cons cs;
  cs.push_back(A >= 0);
  cs.push_back(B >= 0);
  cs.push_back(A < 4);
  cs.push_back(B <= 4);
  cs.push_back(A + B > 0);

  Poly ph1(2, Topol::NNC);
  ph1.add_cons(cs);

  Con c(B > 2);

  Poly ph2 = ph1.split(c);

  Poly knr1(2, Topol::NNC);
  knr1.add_cons(cs);
  knr1.add_con(B > 2);
  knr1.minimize();

  Poly knr2(2, Topol::NNC);
  knr2.add_cons(cs);
  knr2.add_con(B <= 2);
  knr2.minimize();

  return (ph1 == knr1 && ph2 == knr2);
}

bool
test12() {
  Var A(0);
  Var B(1);
  Cons cs;
  cs.push_back(A >= 0);
  cs.push_back(B >= 0);
  cs.push_back(A < 4);
  cs.push_back(B < 4);

  Poly ph1(2, Topol::NNC);
  ph1.add_cons(cs);

  Con c(A + B <= 4);

  Poly ph2 = ph1.split(c);

  Poly knr1(2, Topol::NNC);
  knr1.add_cons(cs);
  knr1.add_con(A + B <= 4);
  knr1.minimize();

  Poly knr2(2, Topol::NNC);
  knr2.add_cons(cs);
  knr2.add_con(A + B > 4);
  knr2.minimize();

  return (ph1 == knr1 && ph2 == knr2);
}

bool
test13() {
  Var A(0);
  Var B(1);
  Cons cs;
  cs.push_back(A >= 0);
  cs.push_back(B >= 0);
  cs.push_back(A < 4);
  cs.push_back(B < 4);

  Poly ph1(2, Topol::NNC);
  ph1.add_cons(cs);

  Con c(A + B <= 0);

  Poly ph2 = ph1.split(c);

  Poly knr1(2, Topol::NNC);
  knr1.add_cons(cs);
  knr1.add_con(A + B <= 0);
  knr1.minimize();

  Poly knr2(2, Topol::NNC);
  knr2.add_cons(cs);
  knr2.add_con(A + B > 0);
  knr2.minimize();

  return (ph1 == knr1 && ph2 == knr2);
}

bool
test14() {
  Var A(0);
  Var B(1);
  Cons cs;
  cs.push_back(A >= 0);
  cs.push_back(B >= 0);
  cs.push_back(A < 4);
  cs.push_back(B < 4);

  Poly ph1(2, Topol::NNC);
  ph1.add_cons(cs);

  Con c(A + B > 4);

  Poly ph2 = ph1.split(c);

  Poly knr1(2, Topol::NNC);
  knr1.add_cons(cs);
  knr1.add_con(A + B > 4);
  knr1.minimize();

  Poly knr2(2, Topol::NNC);
  knr2.add_cons(cs);
  knr2.add_con(A + B <= 4);
  knr2.minimize();

  return (ph1 == knr1 && ph2 == knr2);
}

bool
test15() {
  Var A(0);
  Var B(1);
  Cons cs;
  cs.push_back(A >= 0);
  cs.push_back(B >= 0);
  cs.push_back(A < 4);
  cs.push_back(B < 4);

  Poly ph1(2, Topol::NNC);
  ph1.add_cons(cs);

  Con c(A < 4);

  Poly ph2 = ph1.split(c);

  Poly knr1(2, Topol::NNC);
  knr1.add_cons(cs);
  knr1.add_con(A < 4);
  knr1.minimize();

  Poly knr2(2, Topol::NNC);
  knr2.add_cons(cs);
  knr2.add_con(A >= 4);
  knr2.minimize();

  return (ph1 == knr1 && ph2 == knr2);
}

bool
test16() {
  Var A(0);
  Var B(1);
  Cons cs;
  cs.push_back(A >= 0);
  cs.push_back(A < 4);

  Poly ph1(2, Topol::NNC);
  ph1.add_cons(cs);

  Con c(B <= 4);

  Poly ph2 = ph1.split(c);

  Poly knr1(2, Topol::NNC);
  knr1.add_cons(cs);
  knr1.add_con(B <= 4);

  Poly knr2(2, Topol::NNC);
  knr2.add_cons(cs);
  knr2.add_con(B > 4);

  return (ph1 == knr1 && ph2 == knr2);
}

bool
test17() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Topol::NNC);

  Con c(A + B <= 4);

  Poly ph2 = ph1.split(c);

  Poly knr1(2, Topol::NNC);
  knr1.add_con(A + B <= 4);

  Poly knr2(2, Topol::NNC);
  knr2.add_con(A + B > 4);

  return (ph1 == knr1 && ph2 == knr2);
}

bool
test18() {
  Var A(0);
  Var B(1);
  Cons cs;
  cs.push_back(A >= 0);
  cs.push_back(A <= 3);

  Poly ph1(2, Topol::CLOSED);
  ph1.add_cons(cs);

  Con c(B <= 0);

  Poly ph2 = ph1.split(c);

  Poly knr1(2, Topol::CLOSED);
  knr1.add_cons(cs);
  knr1.add_con(B <= 0);
  knr1.minimize();

  Poly knr2(2, Topol::CLOSED);
  knr2.add_cons(cs);
  knr2.add_con(B >= 0);
  knr2.minimize();

  return (ph1 == knr1 && ph2 == knr2);
}

bool
test19() {
  Var A(0);
  Var B(1);
  Cons cs;
  cs.push_back(A > 0);
  cs.push_back(A < 3);

  Poly ph1(2, Topol::NNC);
  ph1.add_cons(cs);

  Con c(A < 2);

  Poly ph2 = ph1.split(c);

  Poly knr1(2, Topol::NNC);
  knr1.add_cons(cs);
  knr1.add_con(A < 2);
  knr1.minimize();

  Poly knr2(2, Topol::NNC);
  knr2.add_cons(cs);
  knr2.add_con(A >= 2);
  knr2.minimize();

  return (ph1 == knr1 && ph2 == knr2);
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
END_MAIN

