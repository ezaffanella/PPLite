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

  Gens gs1;
  gs1.push_back(point());
  gs1.push_back(ray(B));
  gs1.push_back(ray(-A + B));
  Poly ph1(2, Spec_Elem::EMPTY);
  ph1.add_gens(gs1);

  Gens gs2;
  gs2.push_back(point());
  gs2.push_back(ray(B));
  gs2.push_back(ray(-A - B));
  Poly ph2(2, Spec_Elem::EMPTY);
  ph2.add_gens(gs2);

  print_gens(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");

  ph2.widening_assign(ph1, Widen_Impl::BHRZ03);

  Poly known_result(2);
  known_result.add_con(A <= 0);

  bool ok = (ph2 == known_result);

  print_gens(ph2, "*** after ph2.widening_assign(ph1) ***");

  return ok;
}

bool
test02() {
  Var A(0);
  Var B(1);

  Gens gs1;
  gs1.push_back(point());
  gs1.push_back(ray(-A));
  gs1.push_back(ray(-A - B));
  Poly ph1(2, Spec_Elem::EMPTY);
  ph1.add_gens(gs1);

  Gens gs2;
  gs2.push_back(point());
  gs2.push_back(ray(-A));
  gs2.push_back(ray(A - B));
  Poly ph2(2, Spec_Elem::EMPTY);
  ph2.add_gens(gs2);

  print_gens(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");

  ph2.widening_assign(ph1, Widen_Impl::BHRZ03);

  Poly known_result(2);
  known_result.add_con(B <= 0);

  bool ok = (ph2 == known_result);

  print_gens(ph2, "*** after ph2.widening_assign(ph1) ***");

  return ok;
}

bool
test03() {
  Var A(0);
  Var B(1);

  Gens gs1;
  gs1.push_back(point());
  gs1.push_back(ray(-B));
  gs1.push_back(ray(A - B));
  Poly ph1(2, Spec_Elem::EMPTY);
  ph1.add_gens(gs1);

  Gens gs2;
  gs2.push_back(point());
  gs2.push_back(ray(-B));
  gs2.push_back(ray(A + B));
  Poly ph2(2, Spec_Elem::EMPTY);
  ph2.add_gens(gs2);

  print_gens(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");

  ph2.widening_assign(ph1, Widen_Impl::BHRZ03);

  Poly known_result(2);
  known_result.add_con(A >= 0);

  bool ok = (ph2 == known_result);

  print_gens(ph2, "*** after ph2.widening_assign(ph1) ***");

  return ok;
}

bool
test04() {
  Var A(0);
  Var B(1);

  Gens gs1;
  gs1.push_back(point());
  gs1.push_back(ray(-A));
  gs1.push_back(ray(-A + B));
  Poly ph1(2, Spec_Elem::EMPTY);
  ph1.add_gens(gs1);

  Gens gs2;
  gs2.push_back(point());
  gs2.push_back(ray(-A));
  gs2.push_back(ray(A + B));
  Poly ph2(2, Spec_Elem::EMPTY);
  ph2.add_gens(gs2);

  print_gens(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");

  ph2.widening_assign(ph1, Widen_Impl::BHRZ03);

  Poly known_result(2);
  known_result.add_con(B >= 0);

  bool ok = (ph2 == known_result);

  print_gens(ph2, "*** after ph2.widening_assign(ph1) ***");

  return ok;
}

bool
test05() {
  Var A(0);
  Var B(1);

  Gens gs1;
  gs1.push_back(point());
  gs1.push_back(ray(-B));
  gs1.push_back(ray(-A - B));
  Poly ph1(2, Spec_Elem::EMPTY);
  ph1.add_gens(gs1);

  Gens gs2;
  gs2.push_back(point());
  gs2.push_back(ray(-B));
  gs2.push_back(ray(-A + B));
  Poly ph2(2, Spec_Elem::EMPTY);
  ph2.add_gens(gs2);

  print_gens(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");

  ph2.widening_assign(ph1, Widen_Impl::BHRZ03);

  Poly known_result(2);
  known_result.add_con(A <= 0);

  bool ok = (ph2 == known_result);

  print_gens(ph2, "*** after ph2.widening_assign(ph1) ***");

  return ok;
}

bool
test06() {
  Var A(0);
  Var B(1);

  Gens gs1;
  gs1.push_back(point());
  gs1.push_back(ray(A));
  gs1.push_back(ray(A - B));
  Poly ph1(2, Spec_Elem::EMPTY);
  ph1.add_gens(gs1);

  Gens gs2;
  gs2.push_back(point());
  gs2.push_back(ray(A));
  gs2.push_back(ray(-A - B));
  Poly ph2(2, Spec_Elem::EMPTY);
  ph2.add_gens(gs2);

  print_gens(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");

  ph2.widening_assign(ph1, Widen_Impl::BHRZ03);

  Poly known_result(2);
  known_result.add_con(B <= 0);

  bool ok = (ph2 == known_result);

  print_gens(ph2, "*** after ph2.widening_assign(ph1) ***");

  return ok;
}

bool
test07() {
  Var A(0);
  Var B(1);

  Gens gs1;
  gs1.push_back(point());
  gs1.push_back(ray(B));
  gs1.push_back(ray(A + B));
  Poly ph1(2, Spec_Elem::EMPTY);
  ph1.add_gens(gs1);

  Gens gs2;
  gs2.push_back(point());
  gs2.push_back(ray(B));
  gs2.push_back(ray(A - B));
  Poly ph2(2, Spec_Elem::EMPTY);
  ph2.add_gens(gs2);

  print_gens(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");

  ph2.widening_assign(ph1, Widen_Impl::BHRZ03);

  Poly known_result(2);
  known_result.add_con(A >= 0);

  bool ok = (ph2 == known_result);

  print_gens(ph2, "*** after ph2.widening_assign(ph1) ***");

  return ok;
}

bool
test08() {
  Var A(0);
  Var B(1);

  Poly ph1(Topol::NNC, 2);
  ph1.add_con(B >= 0);
  ph1.add_con(A + B > 0);
  ph1.add_con(A - B < 1);

  Poly ph2(Topol::NNC, 2);
  ph2.add_con(B >= 0);
  ph2.add_con(A > 0);
  ph2.add_con(A < 1);

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  ph1.widening_assign(ph2, Widen_Impl::BHRZ03);

  Poly known_result(Topol::NNC, 2);
  known_result.add_con(B >= 0);

  bool ok = (ph1 == known_result);

  print_cons(ph1, "*** after widening_assign ***");

  return ok;
}

bool
test09() {
  Var A(0);
  Var B(1);
  Var C(2);

  Gens gs1;
  gs1.push_back(point(A));
  gs1.push_back(closure_point());
  gs1.push_back(ray(A));
  gs1.push_back(ray(B));
  gs1.push_back(ray(A + B + 2*C));
  Poly ph1(Topol::NNC, 3, Spec_Elem::EMPTY);
  ph1.add_gens(gs1);

  Gens gs2;
  gs2.push_back(point(A));
  gs2.push_back(closure_point());
  gs2.push_back(ray(A));
  gs2.push_back(ray(B));
  gs2.push_back(ray(A + B + C));
  Poly ph2(Topol::NNC, 3, Spec_Elem::EMPTY);
  ph2.add_gens(gs2);

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  ph1.widening_assign(ph2, Widen_Impl::BHRZ03);

  Poly known_result(Topol::NNC, 3);
  known_result.add_con(A > 0);
  known_result.add_con(B >= 0);
  known_result.add_con(C >= 0);

  bool ok = (ph1 == known_result);

  print_cons(ph1, "*** after widening_assign ***");

  return ok;
}

Gens
aux1_test10(const Integer& half_side) {
  // Square.
  Var A(0);
  Var B(1);
  Gens gs;
  gs.push_back(point(half_side*A + half_side*B));
  gs.push_back(point(half_side*A - half_side*B));
  gs.push_back(point(-half_side*A - half_side*B));
  gs.push_back(point(-half_side*A + half_side*B));
  return gs;
}

Gens
aux2_test10(const Integer& half_diagonal) {
  // Rhombus.
  Var A(0);
  Var B(1);
  Gens gs;
  gs.push_back(point(half_diagonal*A));
  gs.push_back(point(half_diagonal*B));
  gs.push_back(point(-half_diagonal*A));
  gs.push_back(point(-half_diagonal*B));
  return gs;
}

Poly
aux3_test10(unsigned n) {

  Integer half_diagonal = 2;
  for (unsigned i = n / 8; i-- > 0; ) {
    half_diagonal *= 2;
  }
  Integer half_side = half_diagonal;

  Gens gs;
  if (n % 8 < 4) {
    half_side /= 2;
    gs = aux1_test10(half_side);
    Gens gs2 = aux2_test10(half_diagonal);
    Gens::const_iterator gi = gs2.begin();
    for (int i = n % 8; i-- > 0; )
      gs.push_back(*gi++);
  }
  else {
    gs = aux2_test10(half_diagonal);
    Gens gs2 = aux1_test10(half_side);
    Gens::const_iterator gi = gs2.begin();
    for (int i = n % 8 - 4; i-- > 0; )
      gs.push_back(*gi++);
  }
  Poly ph(2, Spec_Elem::EMPTY);
  ph.add_gens(gs);

  return ph;
}

bool
test10() {
  // Chain condition for widenings:
  // for each increasing chain of descriptions p_0, p_1, ..., p_i, ...,
  // the sequence q_0, q_1, ..., q_i, ... defined by q_0 = p_0 and
  // for each i >= 1, q_i = q_{i-1} \nabla p_i is ultimately stationary.

  // Initialization: set q_0.
  Poly q_i_minus_1 = aux3_test10(0);

  for (unsigned i = 1; i <= 100; ++i) {
    print_gens(q_i_minus_1, "*** Result of the previous iteration ***");

    Poly p_i = aux3_test10(i);
    print_gens(p_i, "*** New stuff ***");

    Poly q_i = q_i_minus_1;
    q_i.poly_hull_assign(p_i);
    print_gens(q_i, "*** Poly-hull of previous with new ***");

    q_i.widening_assign(q_i_minus_1, Widen_Impl::BHRZ03);
    print_gens(q_i, "*** Result of widening poly-hull with new ***");

    if (q_i == q_i_minus_1) {
      Poly known_result(2);

      bool ok = (q_i == known_result);

      print_cons(q_i, "*** The cons of the fix point ***");
      print_gens(q_i, "*** The gens of the fix point ***");

      return ok;
    }
    q_i_minus_1 = q_i;
  }
  return false;
}

Gens
aux1_test11() {
  Var A(0);
  Var B(1);
  Gens new_gs;
  new_gs.push_back(point(A));
  new_gs.push_back(point(-B));
  new_gs.push_back(point(-A));
  new_gs.push_back(point(B));
  return new_gs;
}

Gen
aux2_test11(const Gen& p1, const Gen& p2, unsigned magic_number) {
  // Splitting segment.
  auto expr1 = p1.linear_expr();
  auto expr2 = p2.linear_expr();
  const auto& d1 = p1.divisor();
  const auto& d2 = p2.divisor();
  expr1 *= d2;
  expr2 *= d1;
  expr1 += expr2;
  // The divisor for the average is 2 * d1 * d2.
  // by carefully taking a smaller divisor, we obtain a point
  // that won't be redundant in the polyhedron.
  // NOTE: I am not *sure* this dirty kludge of using such
  // a magic number will always succeed.
  return point((magic_number+1)*expr1, magic_number*2*d1*d2);
}

Gens
aux3_test11(const Gens& gs, unsigned magic_number) {
  // Double gens.
  Gens new_gs;
  Gens::const_iterator i = gs.begin();
  Gens::const_iterator gs_end = gs.end();
  while (true) {
    const Gen& g = *i;
    new_gs.push_back(g);
    ++i;
    if (i != gs_end)
      new_gs.push_back(aux2_test11(g, *i, magic_number));
    else {
      // Split the last segment.
      Gens::const_iterator gs_begin = gs.begin();
      new_gs.push_back(aux2_test11(g, *gs_begin, magic_number));
      break;
    }
  }
  return new_gs;
}

Poly
aux4_test11(unsigned n) {

  unsigned needed_vertices = n + 4;

  unsigned magic_number = 1;
  unsigned magic_factor = 4;
  Gens gs = aux1_test11();
  unsigned gs_vertices = 4;

  while (gs_vertices * 2 <= needed_vertices) {
    magic_number *= magic_factor;
    gs = aux3_test11(gs, magic_number);
    gs_vertices *= 2;
  }

  if (gs_vertices < needed_vertices) {
    magic_number *= magic_factor;
    Gens gs2 = aux3_test11(gs, magic_number);
    Gens::const_iterator gs2_i = gs2.begin();
    for ( ; gs_vertices < needed_vertices; ++gs_vertices) {
      // Skip the even indexed vertices of gs2.
      ++gs2_i;
      // Add the odd indexed vertices of gs2.
      gs.push_back(*gs2_i++);
    }
  }

  Poly ph(2, Spec_Elem::EMPTY);
  ph.add_gens(gs);
  return ph;
}

bool
test11() {
  // Chain condition for widenings:
  // for each increasing chain of descriptions p_0, p_1, ..., p_i, ...,
  // the sequence q_0, q_1, ..., q_i, ... defined by q_0 = p_0 and
  // for each i >= 1, q_i = q_{i-1} \nabla p_i is ultimately stationary.

  // Initialization: set q_0.
  Poly q_i_minus_1 = aux4_test11(0);

  for (unsigned i = 1; i <= 100; ++i) {
    print_gens(q_i_minus_1, "*** Result of the previous iteration ***");

    Poly p_i = aux4_test11(i);
    print_gens(p_i, "*** New stuff ***");

    Poly q_i = q_i_minus_1;
    q_i.poly_hull_assign(p_i);
    print_gens(q_i, "*** Poly-hull of previous with new ***");

    q_i.widening_assign(q_i_minus_1, Widen_Impl::BHRZ03);
    print_gens(q_i, "*** Result of widening poly-hull with new ***");

    if (q_i == q_i_minus_1) {
      Poly known_result(2);
      bool ok = (q_i == known_result);
      print_cons(q_i, "*** The cons of the fix point ***");
      print_gens(q_i, "*** The gens of the fix point ***");
      return ok;
    }
    q_i_minus_1 = q_i;
  }
  return false;
}

Gens
aux1_test12() {
  Var A(0);
  Var B(1);
  Var C(2);

  Gens new_gs;
  new_gs.push_back(ray(A + C));
  new_gs.push_back(ray(-B + C));
  new_gs.push_back(ray(-A + C));
  new_gs.push_back(ray(B + C));
  return new_gs;
}

Gen
aux2_test12(const Gen& r1, const Gen& r2, unsigned magic_number) {
  // Splitting facet.
  Var C(2);
  auto expr1 = r1.linear_expr();
  expr1 += r2.linear_expr();
  // NOTE: I am not *sure* this dirty kludge of using such
  // a magic number will always succeed.
  expr1 *= (magic_number + 1);
  expr1 -= C;
  return ray(expr1);
}

Gens
aux3_test12(const Gens& gs, unsigned magic_number) {
  // Double gens.
  Gens new_gs;
  Gens::const_iterator i = gs.begin();
  Gens::const_iterator gs_end = gs.end();
  while (true) {
    const Gen& g = *i;
    new_gs.push_back(g);
    ++i;
    if (i != gs_end)
      new_gs.push_back(aux2_test12(g, *i, magic_number));
    else {
      // Split the last facet.
      Gens::const_iterator gs_begin = gs.begin();
      new_gs.push_back(aux2_test12(g, *gs_begin, magic_number));
      break;
    }
  }
  return new_gs;
}

Poly
aux4_test12(unsigned n) {
  unsigned needed_facets = n + 4;
  unsigned magic_number = 1;
  unsigned magic_factor = 4;
  Gens gs = aux1_test12();
  unsigned gs_facets = 4;

  while (gs_facets * 2 <= needed_facets) {
    magic_number *= magic_factor;
    gs = aux3_test12(gs, magic_number);
    gs_facets *= 2;
  }

  if (gs_facets < needed_facets) {
    magic_number *= magic_factor;
    Gens gs2 = aux3_test12(gs, magic_number);
    Gens::const_iterator gs2_i = gs2.begin();
    for ( ; gs_facets < needed_facets; ++gs_facets) {
      // Skip the even indexed facets of gs2.
      ++gs2_i;
      // Add the odd indexed facets of gs2.
      gs.push_back(*gs2_i++);
    }
  }
  gs.push_back(point());
  Poly ph(3, Spec_Elem::EMPTY);
  ph.add_gens(gs);
  return ph;
}

bool
test12() {
  // Chain condition for widenings:
  // for each increasing chain of descriptions p_0, p_1, ..., p_i, ...,
  // the sequence q_0, q_1, ..., q_i, ... defined by q_0 = p_0 and
  // for each i >= 1, q_i = q_{i-1} \nabla p_i is ultimately stationary.
  //  Var A(0);
  Var B(1);
  Var C(2);

  // Initialization: set q_0.
  Poly q_i_minus_1 = aux4_test12(0);

  for (unsigned i = 1; i <= 100; ++i) {
    print_gens(q_i_minus_1, "*** Result of the previous iteration ***");

    Poly p_i = aux4_test12(i);
    print_gens(p_i, "*** New stuff ***");

    Poly q_i = q_i_minus_1;
    q_i.poly_hull_assign(p_i);
    print_gens(q_i, "*** Poly-hull of previous with new ***");

    q_i.widening_assign(q_i_minus_1, Widen_Impl::BHRZ03);
    print_gens(q_i, "*** Result of widening poly-hull with new ***");

    if (q_i == q_i_minus_1) {
      Poly known_result(3);
      known_result.add_con(-B + C >= 0);
      bool ok = (q_i == known_result);
      print_cons(q_i, "*** The cons of the fix point ***");
      print_gens(q_i, "*** The gens of the fix point ***");
      return ok;
    }
    q_i_minus_1 = q_i;
  }
  return false;
}

bool
test13() {
  Var A(0);
  Var B(1);

  Gens gs1;
  gs1.push_back(point());
  gs1.push_back(point(A + 2*B));
  gs1.push_back(ray(A));
  gs1.push_back(ray(2*A + B));
  Poly ph1(2, Spec_Elem::EMPTY);
  ph1.add_gens(gs1);

  Gens gs2;
  gs2.push_back(point());
  gs2.push_back(point(A + 2*B));
  gs2.push_back(ray(A));
  gs2.push_back(ray(A + B));
  Poly ph2(2, Spec_Elem::EMPTY);
  ph2.add_gens(gs2);

  print_gens(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");

  ph2.widening_assign(ph1, Widen_Impl::BHRZ03);

  Poly known_result(2);
  known_result.add_con(B >= 0);
  known_result.add_con(2*A- B >= 0);

  bool ok = (ph2 == known_result);

  print_gens(ph2, "*** after ph2.widening_assign(ph1) ***");

  return ok;
}

bool
test14() {
  Var A(0);
  Var B(1);

  Gens gs1;
  gs1.push_back(point());
  gs1.push_back(point(A + 3*B));
  gs1.push_back(ray(A));
  gs1.push_back(ray(2*A - B));
  Poly ph1(2, Spec_Elem::EMPTY);
  ph1.add_gens(gs1);

  Gens gs2;
  gs2.push_back(point());
  gs2.push_back(point(A + 3*B));
  gs2.push_back(ray(A + B));
  gs2.push_back(ray(A - B));
  Poly ph2(2, Spec_Elem::EMPTY);
  ph2.add_gens(gs2);

  print_gens(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");

  ph2.widening_assign(ph1, Widen_Impl::BHRZ03);

  Poly known_result(2);
  known_result.add_con(A >= 0);
  known_result.add_con(3*A - B >= 0);

  bool ok = (ph2 == known_result);

  print_gens(ph2, "*** after ph2.widening_assign(ph1) ***");

  return ok;
}

bool
test15() {
  Var A(0);
  Var B(1);

  Gens gs1;
  gs1.push_back(point());
  gs1.push_back(point(6*A - B));
  gs1.push_back(point(6*B));
  gs1.push_back(point(A + 10*B));
  gs1.push_back(ray(A + 2*B));
  Poly ph1(2, Spec_Elem::EMPTY);
  ph1.add_gens(gs1);

  Gens gs2;
  gs2.push_back(point());
  gs2.push_back(point(6*A - B));
  gs2.push_back(point(6*B));
  gs2.push_back(point(A + 10*B));
  gs2.push_back(ray(A + B));
  gs2.push_back(ray(A + 3*B));
  gs2.push_back(point(-4*A + 3*B, 13));
  gs2.push_back(point(-2*A + B, 8));
  gs2.push_back(point(-A + 12*B, 4));

  Poly ph2(2, Spec_Elem::EMPTY);
  ph2.add_gens(gs2);

  print_gens(ph1, "*** ph1 ***");
  print_cons(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");
  print_cons(ph2, "*** ph2 ***");

  ph2.widening_assign(ph1, Widen_Impl::BHRZ03);

  // This is the result of applying H79.
  Gens gs_known_result;
  gs_known_result.push_back(point(-36*A + 6*B, 25));
  gs_known_result.push_back(ray(A + 4*B));
  gs_known_result.push_back(ray(6*A - B));

  Poly known_result(2, Spec_Elem::EMPTY);
  known_result.add_gens(gs_known_result);

  bool ok = (ph2 == known_result);

  print_gens(ph2, "*** after ph2.widening_assign(ph1) ***");
  print_cons(ph2, "*** after ph2.widening_assign(ph1) ***");

  return ok;
}

bool
test16() {
  Var A(0);
  Var B(1);

  Gens gs1;
  gs1.push_back(closure_point());
  gs1.push_back(closure_point(A + B));
  gs1.push_back(point(2*A + B, 2));
  gs1.push_back(ray(A));
  Poly ph1(Topol::NNC, 2, Spec_Elem::EMPTY);
  ph1.add_gens(gs1);

  Gens gs2;
  gs2.push_back(closure_point());
  gs2.push_back(closure_point(A + B));
  gs2.push_back(closure_point(B, 2));
  gs2.push_back(point(2*A + B, 2));
  gs2.push_back(ray(A));
  Poly ph2(Topol::NNC, 2, Spec_Elem::EMPTY);
  ph2.add_gens(gs2);

  print_gens(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");

  ph2.widening_assign(ph1, Widen_Impl::BHRZ03);

  Poly known_result(Topol::NNC, 2);
  known_result.add_con(B > 0);
  known_result.add_con(B < 1);

  bool ok = (ph2 == known_result);

  print_cons(ph2, "*** after  ph2.widening_assign(ph1) ***");

  return ok;
}

bool
test17() {
  Var A(0);
  Var B(1);

  Gens gs;
  gs.push_back(point());
  gs.push_back(point(B));
  gs.push_back(point(A + 2*B));
  gs.push_back(point(A + B));
  Poly ph(Topol::NNC, 2, Spec_Elem::EMPTY);
  ph.add_gens(gs);

  gs.clear();
  gs.push_back(point());
  gs.push_back(point(B));
  gs.push_back(point(A + 2*B));
  gs.push_back(closure_point(A));
  Poly ph1(Topol::NNC, 2, Spec_Elem::EMPTY);
  ph1.add_gens(gs);

  print_cons(ph, "*** ph ***");
  print_cons(ph1, "*** ph1 ***");

  ph1.widening_assign(ph, Widen_Impl::BHRZ03);

  gs.clear();
  gs.push_back(point());
  gs.push_back(point(B));
  gs.push_back(point(A + 2*B));
  gs.push_back(ray(-B));
  Poly known_result(Topol::NNC, 2, Spec_Elem::EMPTY);
  known_result.add_gens(gs);

  bool ok = (ph1 == known_result);

  print_cons(ph1, "*** ph1.widening_assign(ph) ***");

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
END_MAIN
