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
  Con c;
  return c.is_equality()
    && c.space_dim() == 0
    && c.linear_expr().is_equal_to(Linear_Expr())
    && c.inhomo_term() == 0;
}

bool
test02() {
  Con c = Con::zero_dim_false();
  return c.is_equality()
    && c.space_dim() == 0
    && c.linear_expr().is_equal_to(Linear_Expr())
    && c.inhomo_term() == 1;
}

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
END_MAIN

