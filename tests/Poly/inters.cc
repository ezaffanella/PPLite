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
  Var y(1);
  Var z(2);

  Poly icosahedron(3);
  icosahedron.add_con(4*x - 2*y - z >= -14);
  icosahedron.add_con(4*x + 2*y - z >= -2);
  icosahedron.add_con(x + y >= 1);
  icosahedron.add_con(x + y + 2*z >= 5);
  icosahedron.add_con(x >= -1);
  icosahedron.add_con(x + z >= 1);
  icosahedron.add_con(2*x + y -2*z >= -7);
  icosahedron.add_con(x - y + 2*z >= -1);
  icosahedron.add_con(x - y >= -5);
  icosahedron.add_con(2*x - y - 2*z >= -13);
  icosahedron.add_con(-2*x - y + 2*z >= -1);
  icosahedron.add_con(-x + y >= 1);
  icosahedron.add_con(-x + y -2*z >= -7);
  icosahedron.add_con(-4*x + 2*y + z >= 4);
  icosahedron.add_con(-2*x + y + 2*z >= 5);
  icosahedron.add_con(-x >= -1);
  icosahedron.add_con(-x - z >= -5);
  icosahedron.add_con(-4*x - 2*y + z >= -8);
  icosahedron.add_con(-x - y >= -5);
  icosahedron.add_con(-x - y -2*z >= -13);

  Poly column(3);
  column.add_con(y >= 2);
  column.add_con(y <= 4);
  column.add_con(x >= 0);
  column.add_con(x <= 1);

  Poly computed_result = icosahedron;
  computed_result.intersection_assign(column);

  Poly known_result(3);
  known_result.add_con(-4*x - 2*y + z >= -8);
  known_result.add_con(-4*x + 2*y + z >= 4);
  known_result.add_con(-2*x - y + 2*z >= -1);
  known_result.add_con(-2*x + y + 2*z >= 5);
  known_result.add_con(-x - y - 2*z >= -13);
  known_result.add_con(-x - z >= -5);
  known_result.add_con(-x >= -1);
  known_result.add_con(-x + y - 2*z >= -7);
  known_result.add_con(-y >= -4);
  known_result.add_con(y >= 2);
  known_result.add_con(x >= 0);

  bool ok = (computed_result == known_result);

  print_cons(icosahedron, "*** icosahedron ***");
  print_cons(column, "*** column ***");
  print_cons(computed_result, "*** computed_result ***");
  print_cons(known_result, "*** known_result ***");

  return ok;
}

int
aux_test02(const Poly& ph) {
  if (ph.is_empty() || ph.space_dim() == 0)
    return 0;
  const auto& gs = ph.gens();
  return std::count_if(gs.begin(), gs.end(), std::mem_fn(&Gen::is_point));
}

// Intersection of a pyramid with an half-space of variable height.
bool
test02() {
  Var x(0);
  Var y(1);
  Var z(2);

  // This is the height of the pyramid.
  const Integer pyramid_height = 16;

  // We will intersect it with the half-spaces `z <= k' and `z >= k'
  // with k = i*(height/4) for i = -1, 0, 1, ..., 5.
  struct {
    Integer plane_height;
    int num_points_above;
    int num_points_below;
  } ph_nv[]
      = { {-1*(pyramid_height/4), 5, 0},
          { 0*(pyramid_height/4), 5, 4},
          { 1*(pyramid_height/4), 5, 8},
          { 2*(pyramid_height/4), 5, 8},
          { 3*(pyramid_height/4), 5, 8},
          { 4*(pyramid_height/4), 1, 5},
          { 5*(pyramid_height/4), 0, 5}
      };

  Gens gs;
  gs.push_back(point(0*x + 0*y + 0*z));
  gs.push_back(point(2*x + 0*y + 0*z));
  gs.push_back(point(0*x + 2*y + 0*z));
  gs.push_back(point(2*x + 2*y + 0*z));
  gs.push_back(point(x + y + pyramid_height*z));
  Poly pyramid(3, Spec_Elem::EMPTY);
  pyramid.add_gens(gs);

  print_cons(pyramid, "*** pyramid constraints ***");
  print_gens(pyramid, "*** pyramid generators ***");

  bool ok = true;

  for (dim_type i = 0; i <= 6; ++i) {
    // Above.
    Poly hyper_space_above(3);
    hyper_space_above.add_con(z >= ph_nv[i].plane_height);

    Poly computed_result = pyramid;
    computed_result.intersection_assign(hyper_space_above);

    if (ok
        && aux_test02(computed_result) != ph_nv[i].num_points_above)
      ok = false;

    print_cons(hyper_space_above, "*** hyper_space_above ***");
    print_gens(computed_result, "*** computed_result ***");

    // Below.
    Poly hyper_space_below(3);
    hyper_space_below.add_con(z <= ph_nv[i].plane_height);

    computed_result = pyramid;
    computed_result.intersection_assign(hyper_space_below);

    if (ok
        && aux_test02(computed_result) != ph_nv[i].num_points_below)
      ok = false;

    print_cons(hyper_space_below, "*** hyper_space_below ***");
    print_gens(computed_result, "*** computed_result ***");

  }
  return ok;
}

bool
test03() {
  Var x(0);
  Var y(1);

  Poly ph1(2);
  ph1.add_con(x - y >= 0);
  ph1.add_con(x - y <= 1);

  print_cons(ph1, "*** ph1 ***");

  Poly ph2(2);
  ph2.add_con(x >= 0);
  ph2.add_con(y >= 0);
  ph2.add_con(x <= 1);
  ph2.add_con(y <= 1);

  print_cons(ph2, "*** ph2 ***");

  Poly computed_result = ph1;

  computed_result.intersection_assign(ph2);

  Poly known_result(2);
  known_result.add_con(y >= 0);
  known_result.add_con(x - y >= 0);
  known_result.add_con(x <= 1);

  bool ok = (computed_result == known_result);

  print_cons(computed_result, "*** after intersection_assign ***");

  return ok;
}

bool
test04() {
  Var x(0);
  Var y(1);

  Poly ph1(2);
  ph1.add_con(x >= y);
  ph1.add_con(x >= 0);

  Poly ph2(2, Spec_Elem::EMPTY);

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  Poly computed_result1(ph1);
  computed_result1.intersection_assign(ph2);

  Poly computed_result2(ph1.space_dim());
  const auto& ph1_cs = ph1.cons();
  computed_result2.add_cons(ph1_cs.begin(), ph1_cs.end());
  computed_result2.intersection_assign(ph2);

  Poly known_result(2, Spec_Elem::EMPTY);

  bool ok = (computed_result1 == known_result
             && computed_result2 == known_result);

  print_cons(computed_result1, "*** after intersection_assign ***");
  print_cons(computed_result2, "*** after intersection_assign ***");

  return ok;
}

bool
test05() {
  Var x(0);
  Var y(1);

  Poly ph1(2);
  ph1.add_con(x >= y);

  Poly ph2(2, Spec_Elem::EMPTY);
  ph2.add_gen(point());
  ph2.add_gen(line(x));
  ph2.add_gen(ray(y));

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  ph1.intersection_assign(ph2);

  Poly known_result(2);
  known_result.add_con(y >= 0);
  known_result.add_con(x >= y);

  bool ok = (known_result == ph1);

  print_cons(ph1, "*** after intersection_assign ***");

  return ok;
}

bool
test06() {
  Var x(0);
  Var y(1);

  Gens gs1;
  gs1.push_back(point());
  gs1.push_back(point(3*x));
  gs1.push_back(point(3*y));
  gs1.push_back(point(3*x+ 3*y));
  Poly ph1(2, Spec_Elem::EMPTY);
  ph1.add_gens(gs1);

  Gens gs2;
  gs2.push_back(point(x));
  gs2.push_back(point(4*x));
  gs2.push_back(point(x + 3*y));
  gs2.push_back(point(4*x+ 3*y));
  Poly ph2(2, Spec_Elem::EMPTY);
  ph2.add_gens(gs2);

  print_gens(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");

  ph1.intersection_assign(ph2);

  Poly known_result(2, Spec_Elem::EMPTY);
  known_result.add_gen(point(x));
  known_result.add_gen(point(3*x));
  known_result.add_gen(point(x + 3*y));
  known_result.add_gen(point(3*x + 3*y));

  bool ok = (ph1 == known_result);

  print_cons(ph1, "*** after intersection_assign ***");

  return ok;
}

bool
aux_test07(Poly ph1,
           const Poly ph2,
           // Note intentional call-by-value!
           Poly known_result) {

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  ph1.intersection_assign(ph2);

  print_gens(ph1, "*** after intersection_assign ***");

  return ph1 == known_result;
}

bool
test07() {
  Var x(0);
  Var y(1);

  Poly ph1_1(2);
  ph1_1.add_con(x >= 0);
  ph1_1.add_con(y >= 0);
  ph1_1.add_con(x <= 2);
  ph1_1.add_con(y <= 2);
  Poly ph1_2(ph1_1);

  Poly ph2_1(2);
  ph2_1.add_con(x+y <= 0);
  ph2_1.add_con(x+y >= 2);
  Poly ph2_2(ph2_1);
  Poly ph2_3(ph2_1);
  Poly ph2_4(ph2_1);

  bool ok = aux_test07(ph1_1, ph2_1, ph2_1)
    && aux_test07(ph2_2, ph1_2, ph2_2)
    && aux_test07(ph2_3, ph2_4, ph2_3);

  return ok;
}

bool
test08() {
  Var A(0);

  Poly ph1(2, Spec_Elem::EMPTY);
  Poly ph2(2);
  ph2.add_con(A == 0);

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  Poly known_result(ph1);

  ph1.intersection_assign(ph2);

  bool ok = (ph1 == known_result);

  print_cons(ph1, "*** after ph1.intersection_assign(ph2) ***");

  return ok;
}

bool
test09() {
  Poly ph1;
  Poly ph2(0, Spec_Elem::EMPTY);
  ph2.add_gen(point());

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  Poly known_result = ph1;

  ph1.intersection_assign(ph2);

  bool ok = (ph1 == known_result);

  print_cons(ph1, "*** after ph1.intersection_assign(ph2) ***");

  return ok;
}

bool
test10() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::EMPTY);
  Poly ph2(2);
  ph2.add_con(A - B >= 0);

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  Poly known_result(ph1);

  ph1.intersection_assign(ph2);

  bool ok = (ph1 == known_result);

  print_cons(ph1, "*** after ph1.intersection_assign(ph2) ***");

  return ok;
}

bool
test11() {
  Poly ph1;
  Poly ph2;

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  Poly known_result(ph1);

  ph1.intersection_assign(ph2);

  bool ok = (ph1 == known_result);

  print_cons(ph1, "*** after ph1.intersection_assign(ph2) ***");

  return ok;
}

bool
test12() {
  Var A(0);
  Var B(1);

  Poly ph1(2);
  ph1.minimize();
  ph1.add_con(A == B);
  Poly copy_ph1(ph1);

  Poly ph2(2);
  ph2.minimize();
  ph2.add_con(A - B >= 1);
  Poly copy_ph2 = ph2;

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  ph1.intersection_assign(ph2);
  copy_ph1.intersection_assign(copy_ph2);

  bool ok = (ph1 == copy_ph1);

  print_cons(ph1, "*** after intersection_assign ***");
  print_cons(copy_ph1, "*** after intersection_assign ***");

  return ok;
}

bool
test13() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::EMPTY);
  ph1.add_gen(point());
  ph1.minimize();
  ph1.add_gen(line(A + B));
  Poly copy_ph1 = ph1;

  Poly ph2(2, Spec_Elem::EMPTY);
  ph2.add_gen(point());
  ph2.minimize();
  ph2.add_gen(ray(A));
  ph2.add_gen(ray(B));

  Poly copy_ph2 = ph2;

  print_gens(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");

  ph1.intersection_assign(ph2);
  copy_ph1.intersection_assign(copy_ph2);

  bool ok = (ph1 == copy_ph1);

  print_cons(ph1, "*** after intersection_assign ***");
  print_cons(copy_ph1, "*** after intersection_assign ***");

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
END_MAIN
