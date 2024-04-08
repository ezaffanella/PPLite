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

bool
test01() {
  Con c(Var(0) >= 0);
  return check_print(c, "A >= 0");
}

bool
test02() {
  Con c(Var(0) > 1);
  return check_print(c, "A > 1");
}

bool
test03() {
  Con c(Var(0) == -1);
  return check_print(c, "A = -1");
}

bool
test04() {
  Con c(0 <= Var(0));
  return check_print(c, "A >= 0");
}

bool
test05() {
  Con c(1 < Var(0));
  return check_print(c, "A > 1");
}

bool
test06() {
  Con c(-1 == Var(0));
  return check_print(c, "A = -1");
}

bool
test07() {
  Con c(4 <= 2*Var(0));
  return check_print(c, "A >= 2");
}

bool
test08() {
  Con c(1 < 2*Var(0));
  return check_print(c, "2*A > 1");
}

bool
test09() {
  Con c(1 == -3*Var(0));
  return check_print(c, "3*A = -1");
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

