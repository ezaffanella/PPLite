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
  Poly ph_empty(0, Spec_Elem::EMPTY);
  Poly ph_universe(0, Spec_Elem::UNIVERSE);

  // Testing all combinations for 0-dim polyhedra.
  bool ok = true;
  Poly ph;

  // empty, empty
  ph = ph_empty;
  ok &= ph.closed_join_assign_if_exact(ph_empty);
  ok &= (ph == ph_empty);
  print_cons(ph, "*** empty union empty ***");

  // empty, universe
  ph = ph_empty;
  ok &= ph.closed_join_assign_if_exact(ph_universe);
  ok &= (ph == ph_universe);
  print_cons(ph, "*** empty union universe ***");

  // universe, empty
  ph = ph_universe;
  ok &= ph.closed_join_assign_if_exact(ph_empty);
  ok &= (ph == ph_universe);
  print_cons(ph, "*** universe union empty ***");

  // universe, universe
  ph = ph_universe;
  ok &= ph.closed_join_assign_if_exact(ph_universe);
  ok &= (ph == ph_universe);
  print_cons(ph, "*** universe union universe ***");

  return ok;
}

bool
test02() {
  Var x(0);
  Var y(1);

  Poly ph1(2);
  ph1.add_con(x >= -2);
  ph1.add_con(x <= -1);
  ph1.add_con(y >= 0);
  ph1.add_con(y <= 2);

  print_cons(ph1, "*** ph1 ***");

  Poly ph2(2);
  ph2.add_con(x >= 1);
  ph2.add_con(x <= 2);
  ph2.add_con(y >= 0);
  ph2.add_con(y <= 2);

  print_cons(ph2, "*** ph2 ***");

  Poly kr(ph1);

  bool ok = not ph1.closed_join_assign_if_exact(ph2);
  ok &= (ph1 == kr);

  print_cons(ph1, "*** ph1.closed_join_assign_if_exact(ph2) ***");

  return ok;
}

bool
test03() {
  Var x(0);
  Var y(1);

  Poly ph1(2);
  ph1.add_con(x >= -2);
  ph1.add_con(x <= 0);
  ph1.add_con(y >= 0);
  ph1.add_con(y <= 2);

  print_cons(ph1, "*** ph1 ***");

  Poly ph2(2);
  ph2.add_con(x >= 0);
  ph2.add_con(x <= 2);
  ph2.add_con(y >= 0);
  ph2.add_con(y <= 2);

  print_cons(ph2, "*** ph2 ***");

  Poly kr(2);
  kr.add_con(x >= -2);
  kr.add_con(x <= 2);
  kr.add_con(y >= 0);
  kr.add_con(y <= 2);

  bool ok = ph1.closed_join_assign_if_exact(ph2);
  ok &= (ph1 == kr);

  print_cons(ph1, "*** ph1.closed_join_assign_if_exact(ph2) ***");

  return ok;
}

bool
test04() {
  Var x(0);
  Var y(1);

  Poly ph1(2);
  ph1.add_con(x == 0);
  ph1.add_con(y == 0);

  print_cons(ph1, "*** ph1 ***");

  Poly ph2(2);
  ph2.add_con(x >= 0);
  ph2.add_con(x <= 2);
  ph2.add_con(y >= -2);
  ph2.add_con(y <= 2);

  print_cons(ph2, "*** ph2 ***");

  Poly kr(ph2);

  bool ok = ph1.closed_join_assign_if_exact(ph2);
  ok &= (ph1 == kr);

  print_cons(ph1, "*** ph1.closed_join_assign_if_exact(ph2) ***");

  return ok;
}

bool
test05() {
  Var x(0);
  Var y(1);

  Poly ph1(2);
  ph1.add_con(x >= 0);
  ph1.add_con(y == 0);

  print_cons(ph1, "*** ph1 ***");

  Poly ph2(2);
  ph2.add_con(x >= 0);
  ph2.add_con(y >= 2);
  ph2.add_con(y <= 4);

  print_cons(ph2, "*** ph2 ***");

  Poly kr(ph1);

  bool ok = not ph1.closed_join_assign_if_exact(ph2);
  ok &= (ph1 == kr);

  print_cons(ph1, "*** ph1.closed_join_assign_if_exact(ph2) ***");

  return ok;
}

bool
test06() {
  Var x(0);
  Var y(1);

  Poly ph1(2);
  ph1.add_con(x == y);

  print_cons(ph1, "*** ph1 ***");

  Poly ph2(2);
  ph2.add_con(x == 0);

  print_cons(ph2, "*** ph2 ***");

  Poly kr(ph1);

  bool ok = not ph1.closed_join_assign_if_exact(ph2);
  ok &= (ph1 == kr);

  print_cons(ph1, "*** ph1.closed_join_assign_if_exact(ph2) ***");

  return ok;
}

bool
test07() {
  Var x(0);
  Var y(1);

  Poly ph1(2);
  ph1.add_con(x >= y);

  print_cons(ph1, "*** ph1 ***");

  Poly ph2(2);
  ph2.add_con(x >= 0);

  print_cons(ph2, "*** ph2 ***");

  Poly kr(ph1);

  bool ok = not ph1.closed_join_assign_if_exact(ph2);
  ok &= (ph1 == kr);

  print_cons(ph1, "*** ph1.closed_join_assign_if_exact(ph2) ***");

  return ok;
}

bool
test08() {
  Var x(0);
  Var y(1);

  Poly ph1(2);
  ph1.add_con(x >= y);

  print_cons(ph1, "*** ph1 ***");

  Poly ph2(2);
  ph2.add_con(x <= y);

  print_cons(ph2, "*** ph2 ***");

  Poly kr(2);

  bool ok = ph1.closed_join_assign_if_exact(ph2);
  ok &= (ph1 == kr);

  print_cons(ph1, "*** ph1.closed_join_assign_if_exact(ph2) ***");

  return ok;
}

bool
test09() {
  Var w(0);
  Var x(1);
  Var y(2);
  Var z(3);

  Cons cs {
    4*x - 2*y - z + 14 >= 0,
    4*x + 2*y - z + 2 >= 0,
    x + y - 1 >= 0,
    x + y + 2*z - 5 >= 0,
    x + 1 >= 0,
    x + z - 1 >= 0,
    2*x + y -2*z + 7 >= 0,
    x - y + 2*z + 1 >= 0,
    x - y + 5 >= 0,
    2*x - y - 2*z + 13 >= 0,
    -2*x - y + 2*z + 1 >= 0,
    -x + y - 1 >= 0,
    -x + y -2*z + 7 >= 0,
    -4*x + 2*y + z - 4 >= 0,
    -2*x + y + 2*z - 5 >= 0,
    -x + 1 >= 0,
    -x - z + 5 >= 0,
    -4*x - 2*y + z + 8 >= 0,
    -x - y + 5 >= 0,
    -x - y -2*z +13 >= 0
  };

  Poly icosahedron1(4);
  icosahedron1.add_cons(cs);
  icosahedron1.add_con(w >= 0);
  icosahedron1.add_con(w <= 5);

  Poly icosahedron2(4);
  icosahedron2.add_cons(cs);
  icosahedron2.add_con(w >= -7);
  icosahedron2.add_con(w <= 2);

  Poly known_res(4);
  known_res.add_cons(cs);
  known_res.add_con(w >= -7);
  known_res.add_con(w <= 5);

  Poly comp_res = icosahedron1;
  bool ok = comp_res.closed_join_assign_if_exact(icosahedron2);
  ok &= (comp_res == known_res);

  print_cons(icosahedron1, "*** icosahedron1 ***");
  print_cons(icosahedron2, "*** icosahedron2 ***");
  print_cons(comp_res, "*** computed_result ***");
  print_cons(known_res, "*** known result ***");

  return ok;
}

bool
test10() {
  const dim_type sdim = 39;
  std::vector<Var> X;
  // Note: X indices are meant to start from 1, so we have 39+1 entries;
  // first and secnd entry are both for Var(0).
  // X[1] = Var(0), ..., X[i] = Var(i-1), ..., X[39] = Var(38)
  X.push_back(Var(0));
  for (auto i = 0; i != sdim; ++i)
    X.push_back(Var(i));

  Cons cs;
  cs.push_back(X[1] - X[2] - X[3] == 0);
  cs.push_back(Integer("2386907802506363")*X[1] - X[4] == 0);
  cs.push_back(-X[1] >= -80);
  cs.push_back(X[2] - Integer("3152519739159347")*X[14] >= 0);
  cs.push_back(X[6] + X[7] + X[8] + X[9] - X[14] - X[15] == 0);
  cs.push_back(Integer("2386907802506363")*X[6]
               + Integer("2386907802506363")*X[7]
               + Integer("1080863910568919")*X[8]
               + Integer("7746191359077253")*X[9]
               - X[16] == 0);
  cs.push_back(-X[6] + X[10] >= -80);
  cs.push_back(-X[7] + X[11] >= 0);
  cs.push_back(-X[8] + X[12] >= 0);
  cs.push_back(-X[9] + X[13] >= 0);
  cs.push_back(X[22] - X[23] - X[24] - X[25] == 0);
  cs.push_back(Integer("7746191359077253")*X[22] - X[26] == 0);
  cs.push_back(-X[22] >= -500);
  cs.push_back(X[23] - Integer("3152519739159347")*X[36] >= 0);
  cs.push_back(Integer("7746191359077253")*X[28]
               + Integer("7746191359077253")*X[29]
               + Integer("3512807709348987")*X[30]
               + Integer("3332663724254167")*X[31]
               - X[38] == 0);
  cs.push_back(X[28] + X[29] + X[30] + X[31] - X[36] + X[37] + X[39] == 44);
  cs.push_back(-X[28] + X[32] >= -500);
  cs.push_back(-X[29] + X[33] >= 0);
  cs.push_back(-X[30] + X[34] >= 0);
  cs.push_back(-X[31] + X[35] >= 0);
  cs.push_back(Integer("-2661627379775963")*X[10]
               - Integer("2686397177726501")*X[11]
               - Integer("5422333951354077")*X[12]
               - Integer("5469621747441467")*X[13]
               + X[25]
               - Integer("2466846695892189")*X[32]
               - Integer("4996743786567565")*X[33]
               - Integer("5064297780978123")*X[34]
               - Integer("641481471923585")*X[35] >= 0);
  cs.push_back(X[3] - Integer("7854277750134145")*X[22] >= 0);
  cs.push_back(X[15]
               - Integer("7854277750134145")*X[28]
               - Integer("7782220156096217")*X[29]
               - Integer("7782220156096217")*X[30]
               - Integer("7710162562058289")*X[31] >= 0);
  cs.push_back(Integer("-5422333951354077")*X[1] + X[24] >= 0);
  cs.push_back(X[21] >= 2);
  cs.push_back(-X[16] - X[38] >= -300);

  // note: this pushes twice Var(0) >= 0 (no harm)
  for (auto x : X)
    cs.push_back(x >= 0);

  // cs is the common base for ph1, ph2, kr
  Poly ph1(sdim);
  ph1.add_cons(cs);
  Poly ph2 = ph1;
  Poly kr = ph1;

  // minor modifications of ph1 and ph2
  ph1.add_con(X[25] - X[22] <= 5);
  ph2.add_con(X[25] - X[22] >= 1);

  bool ok = ph1.closed_join_assign_if_exact(ph2);
  ok &= (ph1 == kr);
  return ok;
}

bool
test11() {
  const dim_type sdim = 5;
  Cons cs;
  for (dim_type i = 1; i < sdim; ++i) {
    Var x(i);
    cs.push_back(x >= 0);
    cs.push_back(x <= 4);
  }

  Var x(0);

  Poly cube1(sdim);
  cube1.add_cons(cs);
  cube1.add_con(x >= 0);
  cube1.add_con(x <= 4);

  Poly cube2(sdim);
  cube2.add_cons(cs);
  cube2.add_con(x >= 2);
  cube2.add_con(x <= 6);

  bool ok = cube1.closed_join_assign_if_exact(cube2);

  print_gens(cube1, "*** cube1 ***");
  print_gens(cube2, "*** cube2 ***");

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
