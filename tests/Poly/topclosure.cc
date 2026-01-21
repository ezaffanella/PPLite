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
  Var A(0);
  Var B(1);

  Poly ph1(2, Topol::NNC);
  ph1.add_con(A > 1);
  ph1.add_con(A - B > 0);

  Gens gs;
  gs.push_back(point(2*A));
  gs.push_back(closure_point(A + B));
  gs.push_back(ray(-B));
  gs.push_back(ray(A + B));
  Poly ph2(2, Topol::NNC, Spec_Elem::EMPTY);
  ph2.add_gens(gs);

  print_cons(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");

  ph1.topological_closure_assign();
  ph2.topological_closure_assign();

  bool ok = (ph1 == ph2);

  print_cons(ph1, "*** after ph1.topological_closure_assign() ***");
  print_gens(ph2, "*** after ph2.topological_closure_assign() ***");

  return ok;
}

bool
test02() {
  Poly ph1(0, Topol::NNC);
  Poly ph2(2, Topol::NNC, Spec_Elem::EMPTY);

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  Poly known_result1 = ph1;
  Poly known_result2 = ph2;

  ph1.topological_closure_assign();
  ph2.topological_closure_assign();

  bool ok = (ph1 == known_result1 && ph2 == known_result2);

  print_cons(ph1, "*** after ph1.topological_closure_assign() ***");
  print_cons(ph2, "*** after ph2.topological_closure_assign() ***");

  return ok;
}

bool
test03() {
  Var A(0);
  Var B(1);

  Poly ph(2);
  ph.add_con(A - B == 0);
  ph.add_con(A >= 0);

  Poly known_result(ph);

  print_cons(ph, "*** ph ***");

  ph.topological_closure_assign();

  bool ok = (ph == known_result);

  print_cons(ph, "*** after ph.topological_closure_assign() ***");

  return ok;
}

bool
test04() {
  Var A(0);
  Var B(1);

  Poly ph(2, Topol::NNC);
  ph.copy_gens();
  ph.add_con(A > 0);
  ph.add_con(A == B);

  print_cons(ph, "*** ph ***");

  ph.topological_closure_assign();

  Poly known_result(2, Topol::NNC);
  known_result.add_con(A >= 0);
  known_result.add_con(A == B);

  bool ok = (ph == known_result);

  print_cons(ph, "*** after ph.topological_closure_assign() ***");

  return ok;
}

bool
test05() {
  Var A(0);
  Var B(1);

  Poly ph(2, Spec_Elem::EMPTY, Topol::NNC);
  ph.add_gen(point(A));
  ph.minimize();
  ph.add_gen(closure_point());
  ph.add_gen(ray(A));
  ph.add_gen(ray(B));

  print_gens(ph, "*** ph ***");

  ph.topological_closure_assign();

  Poly known_result(2, Topol::NNC);
  known_result.add_con(A >= 0);
  known_result.add_con(B >= 0);

  bool ok = (ph == known_result);

  print_cons(ph, "*** after ph.topological_closure_assign() ***");

  return ok;
}

bool
test06() {
  Var A(0);
  Var B(1);

  Poly ph(2, Topol::NNC);
  ph.add_con(A > 0);
  ph.add_con(B > 0);
  ph.add_con(A < 2);
  ph.add_con(B < 2);
  ph.minimize();

  print_cons(ph, "*** ph  ***");

  ph.topological_closure_assign();

  print_cons(ph, "*** after ph.topological_closure_assign() ***");

  Poly known_result(2, Topol::NNC);
  known_result.add_con(A >= 0);
  known_result.add_con(B >= 0);
  known_result.add_con(A <= 2);
  known_result.add_con(B <= 2);

  bool ok = (ph == known_result);

  return ok;
}

bool
test07() {
  Var A(0);
  Var B(1);

  Poly ph(2, Topol::NNC);
  ph.add_con(A > 0);
  ph.add_con(B > 0);
  ph.minimize();

  print_cons(ph, "*** ph  ***");

  ph.topological_closure_assign();

  print_cons(ph, "*** after ph.topological_closure_assign() ***");

  Poly known_result(2, Topol::NNC);
  known_result.add_con(A >= 0);
  known_result.add_con(B >= 0);

  bool ok = (ph == known_result);

  return ok;
}

bool
test08() {
  Var A(0);
  Var B(1);

  Poly ph(2, Topol::NNC);
  ph.add_con(A > 0);
  ph.add_con(B > 0);
  ph.add_con(A < 2);
  ph.add_con(B < 2);

  print_cons(ph, "*** ph ***");

  ph.topological_closure_assign();

  print_cons(ph, "*** after ph.topological_closure_assign() ***");

  ph.add_gen(closure_point(A + 3*B));
  ph.minimize();

  print_cons(ph, "*** after add_closure_point ***");

  Poly known_result(2, Spec_Elem::EMPTY, Topol::NNC);
  known_result.add_gen(point());
  known_result.add_gen(point(2*A));
  known_result.add_gen(point(2*B));
  known_result.add_gen(point(2*A + 2*B));
  known_result.add_gen(closure_point(A + 3*B));
  known_result.minimize();

  print_cons(known_result, "*** known_result ***");

  bool ok = (ph == known_result);

  return ok;
}

bool
test09() {
  Var A(0);
  Var B(1);

  Poly ph(2, Spec_Elem::EMPTY, Topol::NNC);
  ph.add_gen(point(A + B));
  ph.add_gen(closure_point());
  ph.add_gen(closure_point(2*A));
  ph.add_gen(closure_point(2*B));
  ph.add_gen(closure_point(2*A + 2*B));

  ph.topological_closure_assign();

  print_cons(ph, "*** after ph.topological_closure_assign() ***");

  ph.add_gen(closure_point(A + 3*B));
  ph.minimize();

  print_cons(ph, "*** after add_closure_point ***");

  Poly known_result(2, Spec_Elem::EMPTY, Topol::NNC);
  known_result.add_gen(point());
  known_result.add_gen(point(2*A));
  known_result.add_gen(point(2*B));
  known_result.add_gen(point(2*A + 2*B));
  known_result.add_gen(closure_point(A + 3*B));
  known_result.minimize();

  print_cons(known_result, "*** known_result ***");

  bool ok = (ph == known_result);

  return ok;
}

bool
test10() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Topol::NNC);
  ph1.add_con(A >= 0);
  ph1.add_con(B >= 0);
  ph1.add_con(A + B > 0);

  ph1.topological_closure_assign();
  ph1.set_topology(Topol::CLOSED);

  Poly known_result(2, Spec_Elem::EMPTY);
  known_result.add_gen(point());
  known_result.add_gen(ray(A));
  known_result.add_gen(ray(B));

  bool ok = (ph1 == known_result);
  return ok;
}

bool
test11() {
  Var A(0);
  Var B(1);

  Poly ph1(2);
  ph1.add_con(A >= 0);
  ph1.add_con(B >= 0);
  ph1.set_topology(Topol::NNC);

  ph1.add_con(A + B > 0);

  Poly known_result(2, Spec_Elem::EMPTY, Topol::NNC);
  known_result.add_gen(point(A));
  known_result.add_gen(point(B));
  known_result.add_gen(closure_point());
  known_result.add_gen(ray(A));
  known_result.add_gen(ray(B));

  bool ok = (ph1 == known_result);
  return ok;
}

bool
test12() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Topol::NNC);
  ph1.add_con(A >= 0);
  ph1.add_con(B >= 0);
  ph1.add_con(B < 1);

  ph1.topological_closure_assign();

  return ph1.check_inv();
}

bool
test13() {
  Var A(0);
  Var B(1);

  Poly ph(2, Topol::NNC);
  ph.add_con(A > 0);
  ph.add_con(B >= 0);
  ph.add_con(B <= 1);
  ph.minimize();

  ph.topological_closure_assign();
  ph.add_con(A > 1);
  ph.add_gen(closure_point(2*B));
  ph.minimize();

  Cons cs = ph.copy_cons();
  int ph_strict = std::count_if(cs.begin(), cs.end(),
                                std::mem_fn(&Con::is_strict_inequality));
  return ph_strict == 2;
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
  DO_TEST(test10);
  DO_TEST(test11);
  DO_TEST(test12);
  DO_TEST(test13);
END_MAIN
