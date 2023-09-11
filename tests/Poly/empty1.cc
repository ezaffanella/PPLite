/* PPLite: a lightweight library for convex polyhedra derived from PPL.
   Copyright (C) 2001-2010 Roberto <bagnara@cs.unipr.it>
   Copyright (C) 2010-2018 BUGSENG srl (http://bugseng.com)
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
  Poly ph(0, Spec_Elem::EMPTY);
  return ph.is_necessarily_closed()
    && ph.space_dim() == 0
    && ph.is_empty();
}

bool
test02() {
  Poly ph(5, Spec_Elem::EMPTY);
  return ph.is_necessarily_closed()
    && ph.space_dim() == 5
    && ph.is_empty();
}

bool
test03() {
  Poly ph(2, Spec_Elem::EMPTY, Topol::NNC);
  return !ph.is_necessarily_closed()
    && ph.space_dim() == 2
    && ph.is_empty();
}

bool
test04() {
  Poly ph;
  return ph.is_necessarily_closed()
    && ph.space_dim() == 0
    && !ph.is_empty();
}

bool
test05() {
  Poly ph(2);
  return ph.is_necessarily_closed()
    && ph.space_dim() == 2
    && !ph.is_empty();
}

bool
test06() {
  Poly ph(2);
  ph.add_con(Linear_Expr(0) > 1);
  return ph.is_necessarily_closed()
    && ph.space_dim() == 2
    && ph.is_empty();
}

bool
test07() {
  Poly ph(2);
  ph.add_con(Linear_Expr(0) < 1);
  return ph.is_necessarily_closed()
    && ph.space_dim() == 2
    && !ph.is_empty();
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

