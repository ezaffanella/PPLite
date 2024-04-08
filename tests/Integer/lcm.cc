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

// test Integer lcm functions

bool
test01() {
  // lcm_assign, not prime integers
  Integer i;
  Integer j(15);
  Integer k(12);

  lcm_assign(i, j, k);

  return check_print(i, "60");
}

bool
test02() {
  // lcm_assign, same number
  Integer i;
  Integer j(12);
  Integer k(12);

  lcm_assign(i, j, k);

  return check_print(i, "12");
}

bool
test03() {
  // lcm_assign, 0 and a number
  Integer i;
  Integer j(0);
  Integer k(5);

  lcm_assign(i, j, k);

  return check_print(i, "0");
}

bool
test04() {
  // lcm_assign, 0 and 0
  Integer i;
  Integer j(0);
  Integer k(0);

  lcm_assign(i, j, k);

  return check_print(i, "0");
}

bool
test05() {
  // lcm_assign, prime numbers
  Integer i;
  Integer j(3);
  Integer k(5);

  lcm_assign(i, j, k);

  return check_print(i, "15");
}

bool
test06() {
  // lcm_assign, negatives
  Integer i;
  Integer j(-12);
  Integer k(-16);

  lcm_assign(i, j, k);

  return check_print(i, "48");
}

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
  DO_TEST(test03);
  DO_TEST(test04);
  DO_TEST(test05);
  DO_TEST(test06);
END_MAIN
