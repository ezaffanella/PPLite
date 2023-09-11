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

#include <list>

bool
test01() {
  BBox box1(1);
  box1[0].set_lb(Rational(0));

  BBox box2(1);
  box2[0].set_lb(Rational(1));

  bool ok = box1.contains(box2);
  return ok;
}

bool
test02() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Topol::NNC);
  ph1.add_con(A > 0);
  ph1.add_con(B >= 0);
  ph1.add_con(B < 1);
  ph1.minimize();

  BBox box(2);
  box[0].set_lb(Rational(0));
  box[1].set_lb(Rational(0));
  box[1].set_ub(Rational(1));
  box.maybe_update_volume_info();

  bool ok = (box == ph1.get_bounding_box());
  return ok;
}

bool
test03() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::EMPTY, Topol::NNC);
  ph1.add_gen(point());
  ph1.add_gen(closure_point(A));
  ph1.add_gen(closure_point(A + B));
  ph1.add_gen(closure_point(B));
  ph1.minimize();

  BBox box(2);
  box[0].set_lb(Rational(0));
  box[0].set_ub(Rational(1));
  box[1].set_lb(Rational(0));
  box[1].set_ub(Rational(1));
  box.maybe_update_volume_info();

  bool ok = (box == ph1.get_bounding_box());
  return ok;
}

bool
test04() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Topol::NNC);
  ph1.add_con(A >= 0);
  ph1.add_con(B >= 0);
  ph1.add_con(B <= 1);
  ph1.minimize();

  Poly ph2(2, Topol::NNC);
  ph2.add_con(A >= 1);
  ph2.add_con(A <= 2);
  ph2.add_con(B >= 0);
  ph2.add_con(B <= 1);
  ph2.minimize();

  bool ok = ph1.boxed_contains(ph2);

  return ok;
}

bool
test05() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Topol::NNC);
  ph1.add_con(A > 0);
  ph1.add_con(A <= 2);
  ph1.add_con(B > 0);
  ph1.add_con(B <= 2);
  ph1.minimize();

  Poly ph2(2, Topol::NNC);
  ph2.add_con(A >= 0);
  ph2.add_con(A <= 1);
  ph2.add_con(B >= 0);
  ph2.add_con(B <= 1);
  ph2.minimize();

  bool ok = !ph1.boxed_contains(ph2);

  return ok;
}

bool
test06() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Topol::NNC);
  ph1.add_con(A > 0);
  ph1.add_con(B > 0);
  ph1.add_con(B < 2);
  ph1.minimize();

  Poly ph2(2, Topol::NNC);
  ph2.add_con(A >= 0);
  ph2.add_con(B > 0);
  ph2.add_con(B < 2);
  ph2.minimize();

  BBox box1 = ph1.get_bounding_box();
  BBox box2 = ph2.get_bounding_box();

  bool ok = box1.contains(box2)
    && !ph1.boxed_contains(ph2);

  return ok;
}

bool
test07() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Topol::NNC);
  ph1.add_con(A > 0);
  ph1.add_con(A < 2);
  ph1.add_con(B > 0);
  ph1.add_con(B < 2);
  ph1.minimize();

  Poly ph2(2, Topol::NNC);
  ph2.add_con(A > 0);
  ph2.add_con(A < 1);
  ph2.add_con(B > 0);
  ph2.add_con(B < 1);
  ph2.minimize();

  BBox box1 = ph1.get_bounding_box();
  BBox box2 = ph2.get_bounding_box();

  bool ok = box1.contains(box2)
    && ph1.boxed_contains(ph2);

  return ok;
}

bool
test08() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::EMPTY, Topol::NNC);
  ph1.add_gen(point(3*A));
  ph1.add_gen(point(-3*A));
  ph1.add_gen(point(3*B));
  ph1.add_gen(point(-3*B));
  ph1.minimize();

  Poly ph2(2, Spec_Elem::EMPTY, Topol::NNC);
  ph2.add_gen(point(A));
  ph2.add_gen(point(-A));
  ph2.add_gen(point(B));
  ph2.add_gen(point(-B));
  ph2.minimize();

  BBox box1 = ph1.get_bounding_box();
  BBox box2 = ph2.get_bounding_box();

  bool ok = box1.contains(box2)
    && ph1.boxed_contains(ph2);

  return ok;
}

bool
test09() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::EMPTY, Topol::NNC);
  ph1.add_gen(point(A));
  ph1.add_gen(point(-A));
  ph1.add_gen(point(2*B));
  ph1.add_gen(point(-2*B));
  ph1.minimize();

  Poly ph2(2, Spec_Elem::EMPTY, Topol::NNC);
  ph2.add_gen(point(A));
  ph2.add_gen(point(-A));
  ph2.add_gen(point(B));
  ph2.add_gen(point(-B));
  ph2.minimize();

  BBox box1 = ph1.get_bounding_box();
  BBox box2 = ph2.get_bounding_box();

  bool ok = box1.contains(box2)
    && ph1.boxed_contains(ph2);

  return ok;
}

bool
test10() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::EMPTY, Topol::NNC);
  ph1.add_gen(point(2*A));
  ph1.add_gen(point(2*B));
  ph1.add_gen(ray(A + B));
  ph1.minimize();

  Poly ph2(2, Spec_Elem::EMPTY, Topol::NNC);
  ph2.add_gen(point(2*A + B));
  ph2.add_gen(point(A + 2*B));
  ph2.add_gen(ray(A + B));
  ph2.minimize();

  BBox box1 = ph1.get_bounding_box();
  BBox box2 = ph2.get_bounding_box();

  bool ok = box1.contains(box2)
    && ph1.boxed_contains(ph2);

  return ok;
}

bool
test11() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::EMPTY, Topol::NNC);
  ph1.add_gen(point(B));
  ph1.add_gen(point(A + 2*B));
  ph1.add_gen(point(A));
  ph1.add_gen(ray(A + B));
  ph1.minimize();

  Poly ph2(2, Spec_Elem::EMPTY, Topol::NNC);
  ph2.add_gen(point(A));
  ph2.add_gen(point(A +2*B));
  ph2.add_gen(ray(A));
  ph2.minimize();

  BBox box1 = ph1.get_bounding_box();
  BBox box2 = ph2.get_bounding_box();

  bool ok = box1.contains(box2)
    && !ph1.boxed_contains(ph2);

  return ok;
}

bool
test12() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::EMPTY, Topol::NNC);
  ph1.add_gen(point(2*A));
  ph1.add_gen(point(2*B));
  ph1.add_gen(closure_point(-2*B));
  ph1.add_gen(closure_point(-2*A));
  ph1.minimize();

  Poly ph2(2, Spec_Elem::EMPTY, Topol::NNC);
  ph2.add_gen(point(A));
  ph2.add_gen(point(A - B));
  ph2.add_gen(point(- A - B));
  ph2.minimize();

  BBox box1 = ph1.get_bounding_box();
  BBox box2 = ph2.get_bounding_box();

  bool ok = box1.contains(box2)
    && !ph1.boxed_contains(ph2);

  return ok;
}

bool
test13() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::EMPTY, Topol::NNC);
  ph1.add_gen(point(2*A));
  ph1.add_gen(point(2*B));
  ph1.add_gen(closure_point(-2*B));
  ph1.add_gen(closure_point(-2*A));
  ph1.minimize();

  Poly ph2(2, Spec_Elem::EMPTY, Topol::NNC);
  ph2.add_gen(point(A));
  ph2.add_gen(point(A - B));
  ph2.add_gen(closure_point(- A - B));
  ph2.minimize();

  BBox box1 = ph1.get_bounding_box();
  BBox box2 = ph2.get_bounding_box();

  bool ok = box1.contains(box2)
    && ph1.boxed_contains(ph2);

  return ok;
}

bool
test14() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::EMPTY, Topol::NNC);
  ph1.add_gen(point(2*A));
  ph1.add_gen(point(2*B));
  ph1.add_gen(closure_point(-2*B));
  ph1.add_gen(closure_point(-2*A));
  ph1.minimize();

  Poly ph2(2, Spec_Elem::EMPTY, Topol::NNC);
  ph2.add_gen(point(A + B));
  ph2.add_gen(point(A - B));
  ph2.add_gen(point(- A + B));
  ph2.add_gen(closure_point(- A - B));
  ph2.minimize();

  BBox box1 = ph1.get_bounding_box();
  BBox box2 = ph2.get_bounding_box();

  bool ok = box1.contains(box2)
    && ph1.boxed_contains(ph2);

  return ok;
}

bool
test15() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::EMPTY, Topol::NNC);
  ph1.add_gen(point(10*A));
  ph1.add_gen(point(10*B));
  ph1.add_gen(point(-10*A));
  ph1.add_gen(point(-10*B));
  ph1.minimize();

  Poly ph2(2, Topol::NNC);
  ph2.add_con(7*A - B >= 0);
  ph2.add_con(A - 15*B <= 15);
  ph2.add_con(2*A + B <= 5);
  ph2.minimize();

  BBox box1 = ph1.get_bounding_box();
  BBox box2 = ph2.get_bounding_box();

  bool ok = box1.contains(box2)
    && ph1.boxed_contains(ph2);

  return ok;
}

bool
test16() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::EMPTY, Topol::NNC);
  ph1.add_gen(point(3*A + 3*B));
  ph1.add_gen(point(-3*A + 3*B));
  ph1.add_gen(point(-3*A - 3*B));
  ph1.add_gen(point(3*A - 3*B));
  ph1.minimize();

  Poly ph2(2, Spec_Elem::EMPTY, Topol::NNC);
  ph2.add_gen(point());
  ph2.add_gen(point(5*A));
  ph2.add_gen(point(10*B));
  ph2.add_gen(point(5*A + 10*B));
  ph2.minimize();

  BBox box1 = ph1.get_bounding_box();
  BBox box2 = ph2.get_bounding_box();

  bool ok = !box1.is_disjoint_from(box2);
  if (!ok) return false;
  ok = !box2.is_disjoint_from(box1);
  if (!ok) return false;

  box1.glb_assign(box2);

  ph1.intersection_assign(ph2);
  BBox expected = ph1.get_bounding_box();
  ok = (box1 == expected);

  return ok;
}

bool
test17() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::EMPTY, Topol::NNC);
  ph1.add_gen(point());
  ph1.add_gen(point(2*A));
  ph1.add_gen(point(2*B));
  ph1.minimize();

  Poly ph2(2, Spec_Elem::EMPTY, Topol::NNC);
  ph2.add_gen(point(A + 3*B));
  ph2.add_gen(point(2*A + 3*B));
  ph2.add_gen(ray(B));
  ph2.minimize();

  BBox box1 = ph1.get_bounding_box();
  BBox box2 = ph2.get_bounding_box();

  bool ok = box1.is_disjoint_from(box2);
  if (!ok) return false;
  ok = box2.is_disjoint_from(box1);
  if (!ok) return false;

  box1.glb_assign(box2);
  ok = box1.is_empty();

  return ok;
}

bool
test18() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Topol::NNC);

  Poly ph2(2, Spec_Elem::EMPTY, Topol::NNC);
  ph2.add_gen(point(A + 3*B));
  ph2.add_gen(point(2*A + 3*B));
  ph2.add_gen(ray(B));
  ph2.minimize();

  BBox box1 = ph1.get_bounding_box();
  BBox box2 = ph2.get_bounding_box();

  bool ok = !box1.is_disjoint_from(box2);
  if (!ok) return false;

  box1.glb_assign(box2);
  ok = (box1 == box2);

  return ok;
}

bool
test19() {
  Var A(0);
  Var B(1);

  Poly ph1(3, Spec_Elem::EMPTY, Topol::NNC);

  BBox box1 = ph1.get_bounding_box();

  Gens gs;
  gs.push_back(point(A));
  gs.push_back(point(2*A));
  gs.push_back(ray(B));

  box1.add_gens(gs);
  ph1.add_gens(gs);
  ph1.minimize();

  bool ok = box1 == ph1.get_bounding_box();

  return ok;
}

bool
test20() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Topol::NNC);

  Poly ph2(2, Spec_Elem::EMPTY, Topol::NNC);
  ph2.add_gen(point(A + B));
  ph2.add_gen(point(A));
  ph2.minimize();

  BBox box1 = ph1.get_bounding_box();

  Gens gs;
  gs.push_back(line(A));
  gs.push_back(point(B));
  gs.push_back(point(2*B));

  box1.add_gens(gs);
  ph1.add_gens(gs);
  ph1.minimize();

  bool ok = box1 == ph1.get_bounding_box();

  return ok;
}

bool
test21() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::EMPTY, Topol::NNC);
  ph1.add_gen(point());
  ph1.add_gen(ray(A + B));
  ph1.add_gen(ray(- A + B));
  ph1.minimize();

  BBox box1 = ph1.get_bounding_box();

  Poly ph2(2, Spec_Elem::EMPTY, Topol::NNC);
  ph2.add_gen(point());
  ph2.add_gen(line(A));
  ph2.add_gen(ray(B));
  ph2.minimize();

  BBox box2 = ph2.get_bounding_box();

  bool ok = (box1 == box2);

  return ok;
}

bool
test22() {
  Var A(0);
  Var B(1);
  Var C(2);

  Poly ph1(3, Spec_Elem::EMPTY, Topol::NNC);
  ph1.add_gen(point());
  ph1.add_gen(line(A + C));
  ph1.add_gen(line(A + B));
  ph1.minimize();

  BBox box1 = ph1.get_bounding_box();

  Poly ph2(3, Topol::NNC);
  ph2.minimize();

  BBox box2 = ph2.get_bounding_box();

  bool ok = (box1 == box2);

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
  DO_TEST(test22);
END_MAIN
