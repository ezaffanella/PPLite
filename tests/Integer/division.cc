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

// test Integer division and remainder functions and operators

bool
test01() {
  // divide Integer by Integer
  Integer i(10);
  Integer j(5);
  Integer k = i / j;
  return check_print(k, "2");
}

bool
test02() {
  // divide Integer by signed integer
  Integer i(16);
  Integer k = i / -2;
  k = k / -2L;
  return check_print(k, "4");
}

bool
test03() {
  // divide Integer by unsigner integer
  Integer i(16);
  Integer k = i / 2U;
  k = k / 2UL;
  return check_print(k, "4");
}

bool
test04() {
  // divide assign by Integer
  Integer i(10);
  Integer j(5);
  i /= j;
  return check_print(i, "2");
}

bool
test05() {
  // divide assign by signed integer
  Integer i(1000);
  i /= 10;
  i /= 10L;
  return check_print(i, "10");
}

bool
test06() {
  // divide assign by unsigned integer
  Integer i(1000);
  i /= 10U;
  i /= 10UL;
  return check_print(i, "10");
}

bool
test07() {
  // going more than LONG_MAX and back
  Integer i;
  Integer j(LONG_MAX);
  i = j * 2;
  i /= LONG_MAX;
  return check_print(i, "2");
}

bool
test08() {
  // remainder assign
  Integer i(44);
  i %= 7;
  return check_print(i, "2");
}

bool
test09() {
  // remainder assign
  Integer i(0);
  i %= 21;
  return check_print(i, "0");
}

bool
test10() {
  // remainder assign
  Integer i(-44);
  i %= 7;
  return check_print(i, "-2");
}

bool
test11() {
  // remainder
  Integer i(44);
  Integer j(7);
  Integer k = i % j;
  return check_print(k, "2");
}

bool
test12() {
  // remainder
  Integer i(0);
  Integer j(21);
  Integer k = i % j;
  return check_print(k, "0");
}

bool
test13() {
  // remainder
  Integer i(-44);
  Integer j(7);
  Integer k = i % j;
  return check_print(k, "-2");
}

bool
test14() {
  // divide and assign for 2^exp
  Integer i(1);
  i >>= 1U;
  i >>= 2UL;
  return check_print(i, "0");
}

bool
test15() {
  // divide and assign for 2^exp
  Integer i(8);
  i >>= 1U;
  i >>= 2UL;
  return check_print(i, "1");
}

bool
test16() {
  // divide for 2^exp
  Integer i(1);
  Integer j;
  j = i >> 1U;
  j = j >> 2UL;
  return check_print(j, "0");
}

bool
test17() {
  // divide for 2^exp
  Integer i(8);
  Integer j;
  j = i >> 1U;
  j = j >> 2UL;
  return check_print(j, "1");
}

bool
test18() {
  Integer i;
  Integer j(-12);
  Integer k(-4);

  exact_div_assign(i, j, k);

  return check_print(i, "3");
}

bool
test19() {
//  - exact_div_assign
  Integer i;
  Integer j(-12);
  Integer k(1);

  exact_div_assign(i, j, k);

  return check_print(i, "-12");
}

bool
test20() {
//  - exact_div_assign
  Integer i;
  Integer j(12);
  Integer k(12);

  exact_div_assign(i, j, k);

  return check_print(i, "1");
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
