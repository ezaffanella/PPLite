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

bool
test01() {
  Gen g;
  return g.is_point()
    && g.space_dim() == 0
    && g.linear_expr().is_equal_to(Linear_Expr())
    && g.divisor() == 1;
}

bool
test02() {
  Gen g = point();
  return g.is_point()
    && g.space_dim() == 0
    && g.linear_expr().is_equal_to(Linear_Expr())
    && g.divisor() == 1;
}

bool
test03() {
  Gen g = closure_point();
  return g.is_closure_point()
    && g.space_dim() == 0
    && g.linear_expr().is_equal_to(Linear_Expr())
    && g.divisor() == 1;
}

bool
test04() {
  Gen g = point(Var(3));
  return g.is_point()
    && g.space_dim() == 4
    && g.linear_expr().is_equal_to(Var(3))
    && g.divisor() == 1;
}

bool
test05() {
  Gen g = point(Var(3), 4);
  return g.is_point()
    && g.space_dim() == 4
    && g.linear_expr().is_equal_to(Var(3))
    && g.divisor() == 4;
}

bool
test06() {
  Gen g = point(2*Var(3), 4);
  return g.is_point()
    && g.space_dim() == 4
    && g.linear_expr().is_equal_to(Var(3))
    && g.divisor() == 2;
}

bool
test07() {
  Gen g = ray(-Var(0) + 2 * Var(1));
  return g.is_ray()
    && g.space_dim() == 2
    && g.linear_expr().is_equal_to(-Var(0) + 2 * Var(1));
}

bool
test08() {
  Gen g = line(2 * Var(10));
  return g.is_line()
    && g.space_dim() == 11
    && g.linear_expr().is_equal_to(Var(10));
}

bool
test09() {
  Gen g = line(-2 * Var(5) + 4 * Var(10));
  return g.is_line()
    && g.space_dim() == 11
    && g.linear_expr().is_equal_to(Var(5) - 2 * Var(10));
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

