/* PPLite: a lightweight library for convex polyhedra derived from PPL.
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

// test Integer constructors and assignment operator

bool
test01() {
  // *default* constructor + bool
  Integer i;
  return check_print(i, "0");
}

bool
test02() {
  // static zero member
  return check_print(Integer::zero(), "0");
}

bool
test03() {
  // static one member
  return check_print(Integer::one(), "1");
}

bool
test04() {
  // unsigned integer constructor
  Integer i(12345U);
  Integer j(12345UL);
  bool t1 = check_print(i, "12345");
  bool t2 = check_print(j, "12345");
  return t1 && t2;
}

bool
test05() {
  // signed integer constructor
  Integer i(-1984);
  Integer j(-1984L);
  bool t1 = check_print(i, "-1984");
  bool t2 = check_print(j, "-1984");
  return t1 && t2;
}

bool
test06() {
  // copy constructor
  Integer j(8);
  Integer i(j);
  return check_print(i, "8");
}

bool
test07() {
  // move constructor
  Integer i(Integer(-818));
  return check_print(i, "-818");
}

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
  DO_TEST(test03);
  DO_TEST(test04);
  DO_TEST(test05);
  DO_TEST(test06);
  DO_TEST(test07);
END_MAIN
