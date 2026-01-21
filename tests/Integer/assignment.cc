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

// test Integer assignment operator

bool
test01() {
  // unsigned integer assignment
  Integer i;
  Integer j;
  Integer k;
  i = 12345U;
  j = 12345UL;
  bool t1 = check_print(i, "12345");
  bool t2 = check_print(j, "12345");
  return t1 && t2;
}

bool
test02() {
  // signed integer assignment
  Integer i;
  Integer j;
  Integer k;
  i = -1984;
  j = -1984L;
  bool t1 = check_print(i, "-1984");
  bool t2 = check_print(j, "-1984");
  return t1 && t2;
}

bool
test03() {
  // assignment from another Integer
  Integer j(8);
  Integer i;
  i = j;
  return check_print(i, "8");
}

bool
test04() {
  // move assignment
  Integer i;
  i = Integer(-181);
  return check_print(i, "-181");
}

bool
test05() {
  // testing the swap method
  Integer i(21);
  Integer j(8);
  Integer k;

  swap(i, j);

  bool t1 = check_print(i, "8");
  bool t2 = check_print(j, "21");
  
  swap(j, k);

  bool t3 = check_print(j, "0");
  bool t4 = check_print(k, "21");
  
  return t1 && t2 && t3 && t4;
}


BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
  DO_TEST(test03);
  DO_TEST(test04);
  DO_TEST(test05);
END_MAIN
