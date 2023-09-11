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

bool
test01() {
  Gen g;
  return check_print(g, "p(0)");
}

bool
test02() {
  Var A(0);
  Var C(2);
  Gen g = ray(-2*A + 6*C);
  return check_print(g, "r(-A + 3*C)");
}

bool
test03() {
  Var A(0);
  Var C(2);
  Gen g = line(-2*A + 6*C);
  return check_print(g, "l(A - 3*C)");
}

bool
test04() {
  Var A(0);
  Var C(2);
  Gen g = closure_point(6*A + 6*C, 3);
  return check_print(g, "c(2*A + 2*C)");
}

bool
test05() {
  Var A(0);
  Var C(2);
  Gen g = closure_point(6*A + 6*C, 5);
  return check_print(g, "c((6*A + 6*C)/5)");
}

bool
test06() {
  Var A(0);
  Var C(2);
  Gen g = point(2*A + 6*C, -3);
  return check_print(g, "p((-2*A - 6*C)/3)");
}

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
  DO_TEST(test03);
  DO_TEST(test04);
  DO_TEST(test05);
  DO_TEST(test06);
END_MAIN

