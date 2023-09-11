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

// Testing conversion with combinatorial input.

bool
test01() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::UNIVERSE, Topol::NNC);
  ph1.add_con(A >= 0);
  ph1.add_con(B >= 0);
  ph1.add_con(A <= 2);
  ph1.add_con(B <= 2);
  Index_Set ns1;
  ns1.set(0);
  ns1.set(1);
  ph1.impl().cs_pending.ns_rows.push_back(ns1);

  Poly known_result(2, Spec_Elem::EMPTY, Topol::NNC);
  known_result.add_gen(point(2*B));
  known_result.add_gen(point(2*A));
  known_result.add_gen(point(2*A + 2*B));
  known_result.add_gen(closure_point());

  bool ok = (ph1 == known_result);
  return ok;
}

// Testing redundant indices removal:
// {1,4} becomes {1,3}
bool
test02() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::UNIVERSE, Topol::NNC);
  ph1.add_con(A >= 0);
  ph1.add_con(A <= 2);
  ph1.add_con(A <= 10);
  ph1.add_con(B >= 0);
  ph1.add_con(B <= 2);
  Index_Set ns1;
  ns1.set(0);
  ns1.set(3);
  ph1.impl().cs_pending.ns_rows.push_back(ns1);

  Poly known_result(2, Spec_Elem::EMPTY, Topol::NNC);
  known_result.add_gen(point(2*B));
  known_result.add_gen(point(2*A));
  known_result.add_gen(point(2*A + 2*B));
  known_result.add_gen(closure_point());

  bool ok = (ph1 == known_result);
  return ok;
}

// Adding 2 ns cons.
bool
test03() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::UNIVERSE, Topol::NNC);
  ph1.add_con(A >= 0);
  ph1.add_con(B >= 0);
  ph1.add_con(B <= 2);
  ph1.add_con(A <= 2);
  Index_Set ns1;
  ns1.set(0);
  ns1.set(1);
  ph1.impl().cs_pending.ns_rows.push_back(ns1);
  Index_Set ns2;
  ns2.set(0);
  ns2.set(2);
  ph1.impl().cs_pending.ns_rows.push_back(ns2);

  Poly known_result(2, Spec_Elem::EMPTY, Topol::NNC);
  known_result.add_gen(point(B));
  known_result.add_gen(point(2*A));
  known_result.add_gen(point(2*A + 2*B));
  known_result.add_gen(closure_point());
  known_result.add_gen(closure_point(2*B));

  bool ok = (ph1 == known_result);
  return ok;
}

// Showed bug in simplify
bool
test04() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::UNIVERSE, Topol::NNC);
  ph1.add_con(A >= 0);
  ph1.add_con(B >= 0);
  ph1.add_con(A <= 2);
  ph1.add_con(B <= 2);
  Index_Set ns1;
  ns1.set(2);
  ns1.set(3);
  ph1.impl().cs_pending.ns_rows.push_back(ns1);

  Poly known_result(2, Spec_Elem::EMPTY, Topol::NNC);
  known_result.add_gen(point());
  known_result.add_gen(point(2*B));
  known_result.add_gen(point(2*A));
  known_result.add_gen(closure_point(2*A + 2*B));

  bool ok = (ph1 == known_result);
  return ok;
}

// Adding ns made redundant by another ns cons.
bool
test05() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::UNIVERSE, Topol::NNC);
  ph1.add_con(A >= 0);
  ph1.add_con(B >= 0);
  ph1.add_con(A >= 1);
  ph1.add_con(B >= 1);
  ph1.add_con(A <= 2);
  ph1.add_con(B <= 2);
  Index_Set ns1;
  ns1.set(0);
  ns1.set(1);
  ph1.impl().cs_pending.ns_rows.push_back(ns1);
  Index_Set ns2;
  ns2.set(2);
  ns2.set(3);
  ph1.impl().cs_pending.ns_rows.push_back(ns2);

  Poly known_result(2, Spec_Elem::EMPTY, Topol::NNC);
  known_result.add_gen(point(A + 2*B));
  known_result.add_gen(point(2*A + B));
  known_result.add_gen(point(2*A + 2*B));
  known_result.add_gen(closure_point(A + B));

  bool ok = (ph1 == known_result);
  return ok;
}

// New ns cons create ns gens and
// new ns gens create ns cons
bool
test06() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::UNIVERSE, Topol::NNC);
  ph1.add_con(A >= 0);
  ph1.add_con(B >= 0);
  ph1.add_con(A <= 2);
  ph1.add_con(B <= 2);
  Index_Set ns1;
  ns1.set(0);
  ns1.set(1);
  ph1.impl().cs_pending.ns_rows.push_back(ns1);
  Index_Set ns2;
  ns2.set(0);
  ns2.set(3);
  ph1.impl().cs_pending.ns_rows.push_back(ns2);

  Poly known_result(2, Spec_Elem::EMPTY, Topol::NNC);
  known_result.add_gen(point(2*A));
  known_result.add_gen(point(2*A + 2*B));
  known_result.add_gen(closure_point(2*B));
  known_result.add_gen(closure_point());
  Index_Set ns3;
  ns3.set(1);
  ns3.set(2);
  known_result.impl().gs_pending.ns_rows.push_back(ns3);

  bool ok = (ph1 == known_result);
  return ok;
}

// Adding ns made redundant by a skeleton strict.
bool
test07() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::UNIVERSE, Topol::NNC);
  ph1.add_con(A > 0);
  ph1.add_con(B >= 0);
  ph1.add_con(A <= 2);
  ph1.add_con(B <= 2);
  Index_Set ns1;
  ns1.set(0);
  ns1.set(1);
  ph1.impl().cs_pending.ns_rows.push_back(ns1);

  Poly known_result(2, Spec_Elem::UNIVERSE, Topol::NNC);
  known_result.add_con(A > 0);
  known_result.add_con(A <= 2);
  known_result.add_con(B >= 0);
  known_result.add_con(B <= 2);

  bool ok = (ph1 == known_result);
  return ok;
}

// 3D test:
// ns1 cuts 0-dim face, ns2 cuts 1-dim face
// ns_g1 fills 0-dim face and becomes sk_point
// ns_g2 fills 1-dim face
bool
test08() {
  Var A(0);
  Var B(1);
  Var C(2);

  Poly ph1(3, Spec_Elem::UNIVERSE, Topol::NNC);
  ph1.add_con(A >= 0);
  ph1.add_con(B >= 0);
  ph1.add_con(C >= 0);
  ph1.add_con(A <= 2);
  ph1.add_con(B <= 2);
  ph1.add_con(C <= 2);
  Index_Set ns1;
  ns1.set(0);
  ns1.set(1);
  ns1.set(2);
  ph1.impl().cs_pending.ns_rows.push_back(ns1);
  Index_Set ns2;
  ns2.set(1);
  ns2.set(3);
  ph1.impl().cs_pending.ns_rows.push_back(ns2);

  Poly known_result(3, Spec_Elem::EMPTY, Topol::NNC);
  known_result.add_gen(point(2*B));
  known_result.add_gen(point(2*C));
  known_result.add_gen(point(2*A + 2*B));
  known_result.add_gen(point(2*B + 2*C));
  known_result.add_gen(closure_point(2*A + 2*B + 2*C));
  known_result.add_gen(closure_point());
  known_result.add_gen(closure_point(2*A));
  known_result.add_gen(closure_point(2*A + 2*C));
  Index_Set ns_g1;
  ns_g1.set(3);
  known_result.impl().gs_pending.ns_rows.push_back(ns_g1);
  Index_Set ns_g2;
  ns_g2.set(4);
  ns_g2.set(5);
  known_result.impl().gs_pending.ns_rows.push_back(ns_g2);

  bool ok = (ph1 == known_result);
  return ok;
}

// Adding singleton
bool
test09() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::UNIVERSE, Topol::NNC);
  ph1.add_con(A >= 0);
  ph1.add_con(B >= 0);
  ph1.add_con(A <= 2);
  ph1.add_con(B <= 2);
  Index_Set ns1;
  ns1.set(0);
  ph1.impl().cs_pending.ns_rows.push_back(ns1);

  Poly known_result(2, Spec_Elem::UNIVERSE, Topol::NNC);
  known_result.add_con(A > 0);
  known_result.add_con(A <= 2);
  known_result.add_con(B >= 0);
  known_result.add_con(B <= 2);

  bool ok = (ph1 == known_result);
  return ok;
}

// Incremental add
bool
test10() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::UNIVERSE, Topol::NNC);
  ph1.add_con(A >= 0);
  ph1.add_con(A <= 4);
  ph1.add_con(B >= 0);
  ph1.add_con(B <= 4);
  ph1.minimize();

  ph1.add_con(A <= 2);
  ph1.add_con(B <= 2);
  Index_Set ns1;
  ns1.set(0);
  ns1.set(1);
  ph1.impl().cs_pending.ns_rows.push_back(ns1);

  Poly known_result(2, Spec_Elem::EMPTY, Topol::NNC);
  known_result.add_gen(point());
  known_result.add_gen(point(2*A));
  known_result.add_gen(point(2*B));
  known_result.add_gen(closure_point(2*A + 2*B));

  bool ok = (ph1 == known_result);
  return ok;
}

bool
test11() {
  Var A(0);
  Var B(1);
  Var C(2);

  Poly ph(3, Spec_Elem::UNIVERSE, Topol::NNC);
  ph.add_con(C == 0);
  ph.add_con(A >= 0);
  ph.add_con(B >= 0);
  ph.add_con(B <= 2);
  ph.add_con(A <= 2);
  Index_Set ns1;
  ns1.set(0);
  ns1.set(1);
  ph.impl().cs_pending.ns_rows.push_back(ns1);
  Index_Set ns2;
  ns2.set(0);
  ns2.set(2);
  ph.impl().cs_pending.ns_rows.push_back(ns2);

  ph.minimize();

  Poly known_result(3, Spec_Elem::EMPTY, Topol::NNC);
  known_result.add_gen(point(B));
  known_result.add_gen(point(B));
  known_result.add_gen(point(2*A));
  known_result.add_gen(point(2*A + 2*B));
  known_result.add_gen(closure_point());
  known_result.add_gen(closure_point(2*B));

  bool ok = (ph == known_result);
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
  DO_TEST(test10);
  DO_TEST(test11);
END_MAIN
