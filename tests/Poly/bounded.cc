/* PPLite: a lightweight library for convex polyhedra derived from PPL.
   Copyright (C) 2001-2010 Roberto <bagnara@cs.unipr.it>
   Copyright (C) 2010-2018 BUGSENG srl (http://bugseng.com)
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
  Var x(0);
  Poly ph(2);
  ph.add_con(x >= 0);
  print_cons(ph.copy_cons());
  return !ph.is_bounded();
}

bool
test02() {
  Var x(0);
  Var y(1);
  Poly ph(2);
  ph.add_con(x >= 2);
  ph.add_con(y >= 2);
  ph.add_con(x <= 4);
  ph.add_con(y <= 4);
  print_cons(ph.copy_cons());
  return ph.is_bounded();
}

bool
test03() {
  Poly ph;
  print_cons(ph.copy_cons());
  return ph.is_bounded();
}

bool
test04() {
  Poly ph;
  ph.add_con(Linear_Expr() >= 3);
  print_cons(ph.copy_cons());
  return ph.is_bounded();
}

bool
test05() {
  Poly ph(4, Spec_Elem::EMPTY);
  print_cons(ph.copy_cons());
  return ph.is_bounded();
}

bool
test06() {
  Var A(0), B(1), C(2);
  Poly ph(3, Topol::NNC);
  ph.add_cons({ A + B < 7, C >= 4 });
  print_gens(ph, "=== ph ===");
  bool ok =
    not ph.is_bounded_expr(false, A) &&          // due to line(A - B)
    not ph.is_bounded_expr(true, A) &&           // due to line(A - B)
    not ph.is_bounded_expr(false, B) &&          // due to line(A - B)
    not ph.is_bounded_expr(true, B) &&           // due to line(A - B)
    not ph.is_bounded_expr(false, C) &&          // due to ray(C)
    ph.is_bounded_expr(true, C) &&
    not ph.is_bounded_expr(false, A + B + C) &&  // due to ray(C)
    not ph.is_bounded_expr(true, A + B + C) &&   // due to ray(-A)
    not ph.is_bounded_expr(false, -A - B + C) && // due to ray(C)
    not ph.is_bounded_expr(true, A + B - C);     // due to ray(C)
  return ok;
}

bool
test07() {
  Poly ph(1, Topol::NNC, Spec_Elem::UNIVERSE);
  Var A(0);
  ph.add_con(A > 1);
  return ph.is_bounded_expr(true, A)
    && not ph.is_bounded_expr(false, A)
    && not ph.is_bounded_expr(true, -A)
    && ph.is_bounded_expr(false, -A);
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
