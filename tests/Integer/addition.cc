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

// test Integer sum and difference operators

bool
test01() {
  // sum Integer and Integer
  Integer i;
  Integer j(1);
  Integer k = i + j;
  return check_print(k, "1");
}

bool
test02() {
  // sum Integer and signed integer
  Integer i;
  Integer k = i + 1;
  k = k + 1L;
  return check_print(k, "2");
}

bool
test03() {
  // sum Integer and unsigner integer
  Integer i;
  Integer k = i + 1U;
  k = k + 1UL;
  return check_print(k, "2");
}

bool
test04() {
  // sum assign Integer
  Integer i;
  Integer j(1);
  i += j;
  return check_print(i, "1");
}

bool
test05() {
  // sum assign signed integer
  Integer i;
  i += 1;
  i += 1L;
  return check_print(i, "2");
}

bool
test06() {
  // sum assign unsigned integer
  Integer i;
  i += 1U;
  i += 1UL;
  return check_print(i, "2");
}

bool
test07() {
  // increment by one (pre and post)
  Integer i;
  i++;
  ++i;
  return check_print(i, "2");
}

bool
test08() {
  // add mul, set i = i + j * k
  Integer i(1);
  Integer j(5);
  Integer k(10);
  add_mul_assign(i, j, k);
  return check_print(i, "51");
}

bool
test09() {
  // add mul, set i = i + j * k
  Integer i(-50);
  Integer j(5);
  Integer k(10);
  add_mul_assign(i, j, k);
  return check_print(i, "0");
}

bool
test10() {
  // add mul, set i = i + j * k
  Integer i(-50);
  Integer j(5);
  Integer k(0);
  add_mul_assign(i, j, k);
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
END_MAIN
