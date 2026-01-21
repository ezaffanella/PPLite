/* PPLite: a lightweight library for convex polyhedra derived from PPL.
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

// test Integer multiplication and division operators

bool
test01() {
  // multiply Integer and Integer
  Integer i(2);
  Integer j(3);
  Integer k = i * j;
  return check_print(k, "6");
}

bool
test02() {
  // multiply Integer and signed integer
  Integer i(2);
  Integer k = i * -3;
  k = k * -3L;
  return check_print(k, "18");
}

bool
test03() {
  // multiply Integer and unsigner integer
  Integer i(2);
  Integer k = i * 3U;
  k = k * 3UL;
  return check_print(k, "18");
}

bool
test04() {
  // divide Integer by Integer
  Integer i(10);
  Integer j(5);
  Integer k = i / j;
  return check_print(k, "2");
}

bool
test05() {
  // divide Integer by signed integer
  Integer i(16);
  Integer k = i / -2;
  k = k / -2L;
  return check_print(k, "4");
}

bool
test06() {
  // divide Integer by unsigner integer
  Integer i(16);
  Integer k = i / 2U;
  k = k / 2UL;
  return check_print(k, "4");
}

bool
test07() {
  // multiply assign by Integer
  Integer i(5);
  Integer j(1);
  i *= j;
  return check_print(i, "5");
}

bool
test08() {
  // multiply assign by signed integer
  Integer i(-5);
  i *= -2;
  i *= -2L;
  return check_print(i, "-20");
}

bool
test09() {
  // multiply assign by unsigned integer
  Integer i(-5);
  i *= 2U;
  i *= 2UL;
  return check_print(i, "-20");
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
END_MAIN
