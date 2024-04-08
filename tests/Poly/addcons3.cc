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

namespace {

bool
test01() {
  Var A(0);
  Var B(1);

  Poly ph1(2);
  ph1.add_con(A - B >= 0);
  ph1.add_con(B >= 0);

  Poly ph2(1);
  ph2.add_con(A == 0);

  Cons cs = ph2.copy_cons();

  ph1.add_cons(cs);

  Poly known_result(2);
  known_result.add_con(A == 0);
  known_result.add_con(B == 0);

  bool ok = (ph1 == known_result);

  return ok;
}

bool
test02() {
  Var A(0);
  Var B(1);

  Poly ph1(2);
  ph1.add_con(A >= 0);
  ph1.add_con(B >= 0);

  Cons cs = ph1.copy_cons();

  Poly ph2(3);
  ph2.add_con(A <= 2);

  ph2.add_cons(cs);

  Poly known_result(3);
  known_result.add_con(A >= 0);
  known_result.add_con(A <= 2);
  known_result.add_con(B >= 0);

  bool ok = (ph2 == known_result);

  return ok;
}

bool
test03() {
  set_default_topology(Topol::NNC);

  Var A(0);
  Var B(1);

  Cons cs;
  cs.push_back(A >= 0);
  cs.push_back(B > 0);

  Poly ph1(2);
  ph1.add_cons(cs);
  ph1.minimize();
  return ph1.check_inv();
}

bool
test04() {
  Poly ph1(1, Topol::NNC);
  ph1.add_con(Var(0) > 0);
  ph1.minimize();

  return ph1.check_inv();
}

} // namespace

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
  DO_TEST(test03);
  DO_TEST(test04);
END_MAIN
