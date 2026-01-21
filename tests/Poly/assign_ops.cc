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

// Intersection assign
bool
test01() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::UNIVERSE, Topol::NNC);
  ph1.add_con(A >= 0);
  ph1.add_con(B >= 0);

  Poly ph2(2, Spec_Elem::UNIVERSE, Topol::NNC);
  ph2.intersection_assign(ph1);

  Poly known_result(2, Spec_Elem::EMPTY, Topol::NNC);
  known_result.add_gen(point());
  known_result.add_gen(ray(A));
  known_result.add_gen(ray(B));

  bool ok = (ph2 == known_result
             && ph1 == ph2);
  return ok;
}

// Add ns in intersection
bool
test02() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::UNIVERSE, Topol::NNC);
  ph1.add_con(A >= 0);
  ph1.add_con(B >= 0);
  ph1.add_con(A <= 4);
  ph1.add_con(B <= 4);

  Poly ph2(2, Spec_Elem::EMPTY, Topol::NNC);
  ph2.add_gen(point(2*B));
  ph2.add_gen(point(2*A));
  ph2.add_gen(point(2*A + 2*B));
  ph2.add_gen(closure_point());

  ph1.intersection_assign(ph2);

  bool ok = (ph1 == ph2);
  return ok;
}

bool
test03() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::UNIVERSE, Topol::NNC);
  ph1.add_con(A >= 0);
  ph1.add_con(B >= 0);
  ph1.add_con(A <= 4);
  ph1.add_con(B <= 4);

  Poly ph2(2, Spec_Elem::EMPTY, Topol::NNC);
  ph2.add_gen(point(2*B));
  ph2.add_gen(point(2*A));
  ph2.add_gen(point(2*A + 2*B));
  ph2.add_gen(closure_point());

  ph2.poly_hull_assign(ph1);

  bool ok = (ph1 == ph2);
  return ok;
}

bool
test04() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::UNIVERSE, Topol::NNC);
  ph1.add_con(A >= 1);
  ph1.add_con(B >= 1);
  ph1.add_con(A <= 4);
  ph1.add_con(B <= 4);

  Poly ph2(2, Spec_Elem::EMPTY, Topol::NNC);
  ph2.add_gen(point());
  ph2.add_gen(point(2*A));
  ph2.add_gen(point(2*B));
  ph2.add_gen(closure_point(2*A + 2*B));

  ph2.intersection_assign(ph1);

  Poly known_result(2, Spec_Elem::EMPTY, Topol::NNC);
  known_result.add_gen(point(A + B));
  known_result.add_gen(point(A + 2*B));
  known_result.add_gen(point(2*A + B));
  known_result.add_gen(closure_point(2*A + 2*B));

  bool ok = (ph2 == known_result
             && ph1.check_inv());
  return ok;
}

bool
test05() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::UNIVERSE, Topol::NNC);
  ph1.add_con(A >= 1);
  ph1.add_con(B >= 1);
  ph1.add_con(A <= 4);
  ph1.add_con(B <= 4);

  Poly ph2(2, Spec_Elem::EMPTY, Topol::NNC);
  ph2.add_gen(point());
  ph2.add_gen(point(2*A));
  ph2.add_gen(point(2*B));
  ph2.add_gen(closure_point(2*A + 2*B));

  ph1.intersection_assign(ph2);

  Poly known_result(2, Spec_Elem::EMPTY, Topol::NNC);
  known_result.add_gen(point(A + B));
  known_result.add_gen(point(A + 2*B));
  known_result.add_gen(point(2*A + B));
  known_result.add_gen(closure_point(2*A + 2*B));

  bool ok = (ph1 == known_result
             && ph2.check_inv());
  return ok;
}

bool
test06() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::EMPTY, Topol::NNC);
  ph1.add_gen(point(4*A + 4*B));
  ph1.add_gen(point(4*A + 2*B));
  ph1.add_gen(point(2*B + 4*B));
  ph1.add_gen(closure_point(2*A + 2*B));

  Poly ph2(2, Spec_Elem::EMPTY, Topol::NNC);
  ph2.add_gen(point());
  ph2.add_gen(point(2*A));
  ph2.add_gen(point(2*B));
  ph2.add_gen(closure_point(2*A + 2*B));

  ph1.intersection_assign(ph2);
  ph2.intersection_assign(ph1);

  bool ok = (ph1.is_empty() && ph2.is_empty());
  return ok;
}

bool
test07() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::EMPTY, Topol::NNC);
  ph1.add_gen(point(4*A + 4*B));
  ph1.add_gen(point(4*A));
  ph1.add_gen(point(4*B));
  ph1.add_gen(closure_point());

  Poly ph2(2, Spec_Elem::EMPTY, Topol::NNC);
  ph2.add_gen(point());
  ph2.add_gen(point(2*A));
  ph2.add_gen(point(2*B));
  ph2.add_gen(closure_point(2*A + 2*B));

  ph1.poly_hull_assign(ph2);
  ph2.poly_hull_assign(ph1);

  Poly known_result(2, Spec_Elem::EMPTY, Topol::NNC);
  known_result.add_gen(point());
  known_result.add_gen(point(4*A));
  known_result.add_gen(point(4*B));
  known_result.add_gen(point(4*A + 4*B));

  bool ok = (ph1 == known_result &&
             ph2 == known_result);
  return ok;
}

bool
test08() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::UNIVERSE, Topol::NNC);
  ph1.add_con(2*A - B >= 0);
  ph1.add_con(A - 2*B <= 0);

  Poly ph2(2, Spec_Elem::EMPTY, Topol::NNC);
  ph2.add_gen(point());
  ph2.add_gen(ray(A));
  ph2.add_gen(ray(B));
  ph2.add_con(A + B > 0);

  ph1.intersection_assign(ph2);

  Poly known_result(2, Spec_Elem::UNIVERSE, Topol::NNC);
  known_result.add_con(2*A - B >= 0);
  known_result.add_con(A - 2*B <= 0);
  known_result.add_con(A + B > 0);

  bool ok = (ph1 == known_result
             && ph2.check_inv());
  return ok;
}

// Recognizing efc from ns input
bool
test09() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::EMPTY, Topol::NNC);
  ph1.add_gen(point(A));
  ph1.add_gen(point(2*A));
  ph1.add_gen(line(B));
  // ph1_efc = { A>=1, A<=2 }

  Poly ph2(2, Spec_Elem::EMPTY, Topol::NNC);
  ph2.add_gen(point());
  ph2.add_gen(point(4*A));
  ph2.add_gen(ray(B));
  // In gen2con efc maximal repr is ensured.
  // ph2_efc = { A>=0, A<=4, B>=0 }

  ph1.intersection_assign(ph2);
  ph1.minimize();
  // ph2_efc is not introduced by conversion
  // (even if it were, would be made redundant by ph1_efc)
  // After con2gen: simplify does not ensure maximal efc repr.
  // ph1_efc = { A>=1, A<=2 }

  Poly known_result(2, Spec_Elem::EMPTY, Topol::NNC);
  known_result.add_gen(point(A));
  known_result.add_gen(point(2*A));
  known_result.add_gen(ray(B));
  known_result.minimize();
  // kr_efc = { A>=1, A<=2, B>=0 }

  bool ok = (ph1 == known_result
             && ph2.check_inv());
  return ok;
}

bool
test11() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::EMPTY, Topol::NNC);
  ph1.add_gen(point(A));
  ph1.add_gen(point(2*A));
  ph1.add_gen(line(B));
  ph1.minimize();

  bool ok = (ph1.check_inv());
  return ok;
}

bool
test12() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::EMPTY, Topol::NNC);
  ph1.add_gen(point());
  ph1.add_gen(point(3*A));
  ph1.add_gen(point(3*B));
  ph1.add_gen(point(3*A + 3*B));

  Poly ph2(2, Spec_Elem::EMPTY, Topol::NNC);
  ph2.add_gen(point(2*A + 2*B));
  ph2.add_gen(ray(-A));
  ph2.add_gen(ray(-B));

  ph1.poly_hull_assign(ph2);

  Poly known_result(2, Spec_Elem::EMPTY, Topol::NNC);
  known_result.add_gen(point(3*A + 3*B));
  known_result.add_gen(ray(-A));
  known_result.add_gen(ray(-B));

  ph1.minimize();
  ph1.ascii_dump(nout);
  known_result.ascii_dump(nout);

  bool ok = (ph1 == known_result
             && ph2.check_inv());
  return ok;
}

bool
test13() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Topol::NNC);
  ph1.add_con(A >= 0);
  ph1.add_con(B >= 0);
  ph1.add_con(A + B > 0);

  Poly ph2(2, Topol::NNC);
  ph2.add_con(A >= 0);
  ph2.add_con(B >= 0);

  ph1.poly_hull_assign(ph2);

  ph1.minimize();

  Poly known_result(2, Spec_Elem::EMPTY, Topol::NNC);
  known_result.add_gen(point());
  known_result.add_gen(ray(A));
  known_result.add_gen(ray(B));

  bool ok = (ph1 == known_result
             && ph2.check_inv());
  return ok;
}

bool
test14() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Topol::NNC);
  ph1.add_con(A >= 0);
  ph1.add_con(B >= 0);
  ph1.add_con(A + B > 0);

  Poly ph2(2, Topol::NNC);
  ph2.add_con(A >= 0);
  ph2.add_con(B >= 0);

  ph1.poly_hull_assign(ph2);

  ph1.set_topology(Topol::CLOSED);

  Poly known_result(2, Spec_Elem::EMPTY);
  known_result.add_gen(point());
  known_result.add_gen(ray(A));
  known_result.add_gen(ray(B));

  bool ok = (ph1 == known_result && ph2.check_inv());
  return ok;
}

bool
test15() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Topol::NNC);
  ph1.add_con(A > 0);
  ph1.add_con(B > 0);
  ph1.add_con(A < 1);
  ph1.add_con(B < 1);

  Poly ph2(2, Topol::NNC);
  ph2.add_con(A >= 0);
  ph2.add_con(B >= 0);
  ph2.add_con(A <= 1);
  ph2.add_con(B <= 1);

  ph1.poly_hull_assign(ph2);
  ph1.set_topology(Topol::CLOSED);

  Poly known_result(2, Spec_Elem::EMPTY);
  known_result.add_gen(point());
  known_result.add_gen(point(A));
  known_result.add_gen(point(B));
  known_result.add_gen(point(A + B));

  bool ok = (ph1 == known_result);
  return ok;
}

bool
test16() {
  Var A(0);
  Var B(1);

  Poly ph1(2);
  ph1.add_con(A >= 0);
  ph1.add_con(B >= 0);

  Poly ph2(2);
  ph2.add_con(A <= 2);
  ph2.add_con(B <= 2);

  Poly ph3(2);
  ph3.add_cons(ph1.copy_cons());

  ph1.set_topology(Topol::NNC);
  ph2.set_topology(Topol::NNC);
  ph3.set_topology(Topol::NNC);

  ph1.intersection_assign(ph2);
  ph2.intersection_assign(ph3);

  ph1.set_topology(Topol::CLOSED);
  ph2.set_topology(Topol::CLOSED);

  Poly known_result(2, Spec_Elem::EMPTY);
  known_result.add_gen(point());
  known_result.add_gen(point(2*A));
  known_result.add_gen(point(2*B));
  known_result.add_gen(point(2*A + 2*B));

  bool ok = (ph1 == known_result && ph2 == known_result);
  return ok;
}

bool
test17() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Topol::NNC);
  ph1.add_con(A >= 0);
  ph1.add_con(B >= 0);

  Poly ph2(2, Topol::NNC);
  ph2.add_con(A > 0);
  ph2.add_con(B > 0);

  ph2.poly_hull_assign(ph1);

  if (ph1.is_topologically_closed())
    ph1.set_topology(Topol::CLOSED);
  if (ph2.is_topologically_closed())
    ph2.set_topology(Topol::CLOSED);

  Poly known_result(2, Spec_Elem::EMPTY);
  known_result.add_gen(point());
  known_result.add_gen(ray(A));
  known_result.add_gen(ray(B));

  bool ok = (ph1 == known_result && ph2 == known_result);
  return ok;
}

bool
test18() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Topol::NNC);
  ph1.add_con(A > 0);
  ph1.add_con(B > 0);
  ph1.add_con(A < 1);
  ph1.add_con(B < 1);

  Poly ph2(2, Topol::NNC);
  ph2.add_con(A >= 0);
  ph2.add_con(B >= 0);
  ph2.add_con(A + B > -2);
  ph2.add_con(A < 2);
  ph2.add_con(B < 2);
  ph2.add_con(A <= 1);
  ph2.add_con(B <= 1);

  ph1.poly_hull_assign(ph2);

  if (ph1.is_topologically_closed())
    ph1.set_topology(Topol::CLOSED);
  if (ph2.is_topologically_closed())
    ph2.set_topology(Topol::CLOSED);

  Poly known_result(2, Spec_Elem::EMPTY);
  known_result.add_gen(point());
  known_result.add_gen(point(A));
  known_result.add_gen(point(B));
  known_result.add_gen(point(A + B));

  bool ok = (ph1 == known_result && ph2 == known_result);
  return ok;
}

bool
test19() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Topol::NNC);
  ph1.add_con(A > 0);
  ph1.add_con(B > 0);
  ph1.add_con(A < 1);
  ph1.add_con(B < 1);
  ph1.minimize();

  Cons c;
  ph1.add_cons(c.begin(), c.end());

  bool ok = (ph1.check_inv());
  return ok;
}

bool
test20() {
  Var A(0);
  Var B(1);

  Poly cph1(2, Spec_Elem::EMPTY);
  cph1.add_gen(point());
  cph1.add_gen(point(A));
  cph1.add_gen(ray(B));
  print_cons(cph1, "** C Ph1 minimized cons **");

  Poly cph2(2, Spec_Elem::EMPTY);
  cph2.add_gen(point());
  cph2.add_gen(point(B));
  cph2.add_gen(ray(A));
  print_cons(cph2, "** C Ph2 minimized cons **");

  Poly ph1(2, Topol::NNC);
  ph1.add_cons(cph1.copy_cons());

  Poly ph2(2, Topol::NNC);
  ph2.add_cons(cph2.copy_cons());

  cph2.intersection_assign(cph1);
  print_cons(cph2, "** C Ph2 cons after intersection - rfc in pending **");
  cph2.minimize();
  print_cons(cph2, "** C Ph2 minimized cons after intersection **");

  ph2.intersection_assign(ph1);
  print_cons(ph2, "** Ph2 cons after intersection - rfc in pending **");
  ph2.minimize();
  print_cons(ph2, "** Ph2 minimized cons after intersection **");

  bool ok = (ph1.check_inv() && ph2.check_inv()
             && cph1.check_inv() && cph2.check_inv());
  return ok;
}

bool
test21() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Topol::NNC, Spec_Elem::EMPTY);
  ph1.add_gen(point());
  ph1.add_gen(ray(A));
  ph1.add_gen(ray(B));
  print_cons(ph1, "** Ph1 minimized cons **");

  Poly ph2(2, Topol::NNC, Spec_Elem::EMPTY);
  ph2.add_gen(point());
  ph2.add_gen(point(B));
  ph2.add_gen(ray(A));
  print_cons(ph2, "** Ph2 minimized cons **");

  ph2.intersection_assign(ph1);
  print_cons(ph2, "** Ph2 cons after intersection - skel pos in pending **");
  ph2.minimize();
  print_cons(ph2, "** Ph2 minimized cons after intersection **");

  bool ok = (ph2.check_inv());
  return ok;
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
/* test10 moved into Poly_Impl */
  DO_TEST(test11);
  DO_TEST(test12);
  DO_TEST(test13);
  DO_TEST(test14);
  DO_TEST(test15);
  DO_TEST(test16);
  DO_TEST(test17);
  DO_TEST(test18);
  DO_TEST(test19);
  DO_TEST(test20);
  DO_TEST(test21);
END_MAIN
