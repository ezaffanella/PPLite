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
  Linear_Expr e = Var(0) + 2*Var(2) + 3*Var(5);
  std::vector<dim_type> tbr = { 1, 5 };
  e.remove_space_dims(tbr.begin(), tbr.end());
  return check_print(e, "A + 2*B");
}

bool
test02() {
  Linear_Expr e = Var(0) + 2*Var(2) + 3*Var(5);
  std::vector<dim_type> tbr = { 1, 5 };
  e.remove_space_dims(tbr.begin() + 1, tbr.end());
  return check_print(e, "A + 2*C");
}

bool
test03() {
  Linear_Expr e = Var(0) + 2*Var(2) + 3*Var(5);
  std::vector<dim_type> tbr = { 1, 2, 3, 4, 5, 6, 7, 8 };
  e.remove_space_dims(tbr.begin(), tbr.end());
  return check_print(e, "A");
}

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
  DO_TEST(test03);
END_MAIN
