/* PPLite: a lightweight library for convex polyhedra derived from PPL.
   Copyright (C) 2001-2010 Roberto <bagnara@cs.unipr.it>
   Copyright (C) 2010-2018 BUGSENG srl (http://bugseng.com)
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
  Var x(0);
  Var y(1);
  Var z(2);
  Var w(3);
  Gens gs;
  gs.push_back(point(0*x + y +0*z + 2*w));
  Poly ph(4, Spec_Elem::EMPTY);
  ph.add_gens(gs);

  print_gens(ph, "*** ph ***");

  // This is the set of the variables that we want to remove.
  Index_Set to_be_removed;
  to_be_removed.set(y.id());
  to_be_removed.set(z.id());
  ph.remove_space_dims(to_be_removed);

  Gens known_result_gs;
  known_result_gs.push_back(point(0*x +2*y));
  Poly known_result(2, Spec_Elem::EMPTY);
  known_result.add_gens(known_result_gs);

  bool ok = (known_result == ph);

  print_gens(ph, "*** ph ***");
  print_gens(known_result, "*** known_result ***");

  return ok;
}

bool
test02() {
  Gens gs;

  // Creating 10 points.
  for (int i = 0; i < 10; i++) {
    Linear_Expr e;
    for (int j = 0; j < 10; j++)
      e += (10*i + j) * Var(j);
    gs.push_back(point(e));
  }

  Poly ph(10, Spec_Elem::EMPTY);
  ph.add_gens(gs);

  print_gens(ph, "*** before ***");

  // This is the set of the variables that we want to remove.
  Index_Set to_be_removed;
  to_be_removed.set(0);
  to_be_removed.set(5);
  to_be_removed.set(3);
  to_be_removed.set(4);
  to_be_removed.set(8);

  ph.remove_space_dims(to_be_removed);

  // Useless, but much clearer.
  gs.clear();

  Var a(0);
  Var b(1);
  Var c(2);
  Var d(3);
  Var e(4);

  Linear_Expr expr01 = (1*a + 2*b + 6*c + 7*d + 9*e);
  Linear_Expr expr10 = 10 * (a + b + c + d + e);

  for (int i = 0; i < 10; i++) {
    Linear_Expr expr = i * expr10 + expr01;
    gs.push_back(point(expr));
  }

  Poly known_result(5, Spec_Elem::EMPTY);
  known_result.add_gens(gs);

  bool ok = (ph == known_result);

  print_gens(ph, "*** after ***");
  print_gens(known_result, "*** known_result ***");

  return ok;
}

bool
test03() {
  Var y(1);
  Var z(2);
  Var w(6);

  // This is the set of the variables that we want to remove.
  Index_Set to_be_removed;
  to_be_removed.set(y.id());
  to_be_removed.set(z.id());
  to_be_removed.set(w.id());

  // A 10-dim space, empty polyhedron.
  Poly ph(10, Spec_Elem::EMPTY);
  ph.remove_space_dims(to_be_removed);

  // A 7-dim space, empty polyhedron.
  Poly known_result(7, Spec_Elem::EMPTY);

  bool ok = (known_result == ph);

  print_cons(ph, "*** ph ***");
  print_cons(known_result, "*** known_result ***");

  return ok;
}

bool
test04() {
  Var x(0);
  Var y(1);
  Var z(2);

  Poly ph1(3);
  ph1.add_con(x >= 3);
  ph1.add_con(x - y >= 0);

  print_cons(ph1, "*** ph1 ***");

  Poly ph2 = ph1;

  print_cons(ph2, "*** ph2 ***");

  // This is the set of the variables that we want to remove.
  Index_Set to_be_removed;
  to_be_removed.set(y.id());
  to_be_removed.set(z.id());
  to_be_removed.set(x.id());

  ph1.remove_space_dims(to_be_removed);
  ph2.remove_higher_space_dims(0);

  bool ok = (ph1 == ph2);

  print_gens(ph1, "*** ph1 after remove_space_dims ***");
  print_gens(ph2, "*** ph2 after remove_higher_space_dims ***");

  return ok;
}

bool
test05() {
  Var A(0);
  Gens gs;
  gs.push_back(point());
  gs.push_back(ray(A));
  Poly ph(1, Spec_Elem::EMPTY);
  ph.add_gens(gs);
  ph.add_con(A >= 2);

  print_cons(ph, "*** ph ***");

  Poly known_result(ph);

  // This is the set of the variables that we want to remove.
  Index_Set to_be_removed;

  ph.remove_space_dims(to_be_removed);

  bool ok = (ph == known_result);

  print_cons(ph,
                    "*** after ph.remove_space_dims(to_be_removed) ***");

  return ok;
}

bool
test06() {
  Var x(0);
  Var y(1);
  Var z(2);

  Poly ph1(4, Topol::NNC);

  ph1.add_con(x - y == 3);
  ph1.add_con(z - x > 4);
  ph1.add_con(y < 6);

  print_cons(ph1, "*** ph1 ***");

  Poly ph2(ph1);

  ph1.remove_higher_space_dims(1);

  // This is the set of the variables that we want to remove.
  Index_Set to_be_removed;
  to_be_removed.set(y.id());
  to_be_removed.set(z.id());
  to_be_removed.set(3);

  ph2.remove_space_dims(to_be_removed);

  bool ok = (ph1 == ph2);

  print_cons(ph1, "*** after remove_higher_space_dims(1) ***");
  print_cons(ph2, "*** after remove_space_dims(to_be_removed) ***");

  return ok;
}

bool
test07() {
  Var x(0);
  Var y(1);
  Var z(2);

  Poly ph(3, Topol::NNC);

  ph.add_con(x >= 1);
  ph.add_con(y >= 1);
  ph.add_con(z >= 1);

  print_gens(ph, "*** ph ***");

  Index_Set to_be_removed;
  to_be_removed.set(x.id());
  to_be_removed.set(z.id());

  ph.remove_space_dims(to_be_removed);

  Poly known_result(1, Topol::NNC);
  known_result.add_con(x >= 1);

  bool ok = (ph == known_result);

  print_cons(ph, "*** ph.remove_space_dims() ***");
  print_cons(known_result, "*** known_result ***");

  return ok;
}

bool
test08() {
  Var A(0);
  Var B(1);

  Poly ph(2, Spec_Elem::EMPTY);
  ph.add_gen(point());
  ph.add_gen(ray(2*A + B));

  ph.minimize();
  print_cons(ph, "*** ph ***");
  print_gens(ph, "*** ph ***");

  Index_Set to_be_removed;
  to_be_removed.set(B.id());

  ph.remove_space_dims(to_be_removed);

  Poly known_result(1);
  known_result.add_con(A >= 0);

  return (ph == known_result);
}

bool
test09() {
  Var A(0);
  Var B(1);
  Var C(2);
  Poly ph1(3, Topol::NNC);
  ph1.add_con(A >= 0);
  ph1.add_con(A <= 1);
  ph1.add_con(B >= 0);
  ph1.add_con(A + B > 0);
  ph1.add_con(3*A + 2*B > 0);
  ph1.add_con(C == 0);

  ph1.remove_space_dims(Index_Set(C.id()));

  Poly kr(2, Topol::NNC);
  kr.add_con(A >= 0);
  kr.add_con(A <= 1);
  kr.add_con(B >= 0);
  kr.add_con(A + B > 0);
  kr.add_con(3*A + 2*B > 0);

  return (ph1 == kr);
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
