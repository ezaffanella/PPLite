/* PPLite: a lightweight library for convex polyhedra derived from PPL.
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

// test Integer sum and difference operators

bool
test01() {
  // subtract Integer and Integer
  Integer i;
  Integer j(1);
  Integer k = i - j;
  return check_print(k, "-1");
}

bool
test02() {
  // subtract Integer and signed integer
  Integer i;
  Integer k = i - 1;
  k = k - 1L;
  return check_print(k, "-2");
}

bool
test03() {
  // subtract Integer and unsigner integer
  Integer i;
  Integer k = i - 1U;
  k = k - 1UL;
  return check_print(k, "-2");
}

bool
test04() {
  // subtract assign Integer
  Integer i;
  Integer j(1);
  i -= j;
  return check_print(i, "-1");
}

bool
test05() {
  // subtract assign signed integer
  Integer i;
  i -= 1;
  i -= 1L;
  return check_print(i, "-2");
}

bool
test06() {
  // subtract assign unsigned integer
  Integer i;
  i -= 1U;
  i -= 1UL;
  return check_print(i, "-2");
}

bool
test07() {
  // decrement by one (pre and post)
  Integer i;
  i--;
  --i;
  return check_print(i, "-2");
}

bool
test08() {
  // going more than LONG_MAX and back
  Integer i;
  Integer j(LONG_MAX);
  i = j + j;
  ++i;
  i = i - j - j;
  return check_print(i, "1");
}

bool
test09() {
  // testing the neg method
  Integer i(21);
  Integer j(-8);
  Integer k;

  bool t1 = check_print(neg(i), "-21");
  bool t2 = check_print(i, "21");
  bool t3 = check_print(neg(j), "8");
  bool t4 = check_print(j, "-8");
  bool t5 = check_print(neg(k), "0");
  bool t6 = check_print(k, "0");

  return t1 && t2 && t3 && t4 && t5 && t6;
}

bool
test10() {
  // testing the unary minus operator
  Integer i(21);
  Integer j(-8);
  Integer k;

  bool t1 = check_print(-i, "-21");
  bool t2 = check_print(i, "21");
  bool t3 = check_print(-j, "8");
  bool t4 = check_print(j, "-8");
  bool t5 = check_print(-k, "0");
  bool t6 = check_print(k, "0");

  return t1 && t2 && t3 && t4 && t5 && t6;
}

bool
test11() {
  // testing the neg_assign method
  Integer i(21);
  Integer j(-8);
  Integer k;

  neg_assign(i);
  neg_assign(j);
  neg_assign(k);

  bool t1 = check_print(i, "-21");
  bool t2 = check_print(j, "8");
  bool t3 = check_print(k, "0");

  return t1 && t2 && t3;
}

bool
test12() {
  // sub mul, set i = i - j * k
  Integer i(1);
  Integer j(5);
  Integer k(10);
  sub_mul_assign(i, j, k);
  return check_print(i, "-49");
}

bool
test13() {
  // sub mul, set i = i - j * k
  Integer i(-50);
  Integer j(5);
  Integer k(10);
  sub_mul_assign(i, j, k);
  return check_print(i, "-100");
}

bool
test14() {
  // sub mul, set i = i - j * k
  Integer i(-50);
  Integer j(5);
  Integer k(0);
  sub_mul_assign(i, j, k);
  return check_print(i, "-50");
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
