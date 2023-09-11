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
  Var x(0);
  Var y(1);

  Gens gs;
  gs.push_back(point());
  gs.push_back(point(x));
  gs.push_back(point(y));
  gs.push_back(point(x + y));

  Poly ph(2, Spec_Elem::EMPTY);
  ph.add_gens(gs);

  print_gens(ph, "*** ph ***");

  ph.add_space_dims(1, true);

  Poly known_result(3, Spec_Elem::EMPTY);
  known_result.add_gen(point());
  known_result.add_gen(point(x));
  known_result.add_gen(point(y));
  known_result.add_gen(point(x + y));

  bool ok = (ph == known_result);

  print_gens(ph, "*** after add_space_dims ***");

  return ok;
}

bool
test02() {
  Poly ph(3, Spec_Elem::EMPTY);

  print_cons(ph, "*** ph ***");

  Poly computed_result1(ph);
  Poly computed_result2(ph);

  computed_result1.add_space_dims(4, true);
  computed_result2.add_space_dims(4);

  Poly known_result(7, Spec_Elem::EMPTY);

  bool ok = (computed_result1 == known_result
             && computed_result2 == known_result);

  print_cons(computed_result1, "*** computed_result1 ***");
  print_cons(computed_result2, "*** computed_result2 ***");

  return ok;
}

bool
test03() {
  Var x(0);
  Var y(1);
  Var z(2);
  Var u(3);
  Var v(4);
  Var w(5);

  Gens gs;
  gs.push_back(point());
  gs.push_back(ray(x + y));

  Poly ph(2, Spec_Elem::EMPTY);
  ph.add_gens(gs);

  print_gens(ph, "*** ph ***");

  ph.minimize();

  ph.add_space_dims(2);

  print_gens(ph, "*** after add_space_dims(2) ***");

  ph.add_space_dims(2);

  Poly known_result(6, Spec_Elem::EMPTY);
  known_result.add_gen(point());
  known_result.add_gen(ray(x + y));
  known_result.add_gen(line(z));
  known_result.add_gen(line(u));
  known_result.add_gen(line(v));
  known_result.add_gen(line(w));

  bool ok = (ph == known_result);

  print_gens(ph, "*** ph ***");

  return ok;
}

bool
test04() {
  Poly ph1;

  print_gens(ph1, "*** ph1 ***");

  ph1.add_space_dims(3, true);

  print_gens(ph1, "*** after add_space_dims(3, true) ***");

  Poly ph2;
  print_cons(ph2, "*** ph2 ***");
  print_gens(ph2, "*** ph2 ***");

  ph2.add_space_dims(3, true);

  bool ok = (ph1 == ph2);

  print_gens(ph2, "*** after add_space_dims(3, true) ***");

  return ok;
}

bool
test05() {
  Var C(2);

  Poly ph(2);

  print_cons(ph, "*** ph ***");

  ph.add_space_dims(1, true);

  Poly known_result(3);
  known_result.add_con(C == 0);

  bool ok = (ph == known_result);

  print_cons(ph, "*** after ph.add_space_dims(1, true) ***");

  return ok;
}

bool
test06() {
  Var A(0);

  Poly ph1(2);
  ph1.add_con(A >= 0);
  ph1.add_con(A <= 2);

  Poly ph2(ph1);

  print_cons(ph1, "*** ph1 ***");
  print_cons(ph2, "*** ph2 ***");

  ph1.add_space_dims(0);
  ph2.add_space_dims(0, true);

  bool ok = (ph1 == ph2);

  print_cons(ph1, "*** after ph1.add_space_dims(0) ***");
  print_cons(ph2, "*** after ph2.add_space_dims(0, true) ***");

  return ok;
}


bool
test07() {
  Var x(0);
  Var y(1);
  Var z(2);
  Var w(3);

  Cons cs;
  cs.push_back(x > 2);
  cs.push_back(y > 2);
  cs.push_back(x < 6);
  cs.push_back(y < 6);

  Poly ph(2, Topol::NNC);
  ph.add_cons(cs);

  ph.minimize();

  print_cons(ph, "*** ph ***");
  print_gens(ph, "*** ph ***");

  ph.add_space_dims(2, true);

  Poly known_result(4, Topol::NNC);
  known_result.add_con(z == 0);
  known_result.add_con(w == 0);
  known_result.add_con(x > 2);
  known_result.add_con(y > 2);
  known_result.add_con(x < 6);
  known_result.add_con(y < 6);

  bool ok = (ph == known_result);

  print_cons(ph, "*** after add_space_dims ***");
  print_gens(ph, "*** after add_space_dims ***");

  return ok;
}

bool
test08() {
  Var x(0);

  Poly ph1(1, Topol::NNC);
  ph1.add_con(x > 1);
  ph1.add_con(x < 5);
  print_cons(ph1, "*** ph1 ***");

  ph1.add_space_dims(2);
  ph1.minimize();

  Gens gs;
  gs.push_back(point(2*x));
  gs.push_back(closure_point(x));
  gs.push_back(closure_point(5*x));
  Poly ph2(1, Topol::NNC, Spec_Elem::EMPTY);
  ph2.add_gens(gs);
  print_gens(ph2, "*** ph2 ***");

  ph2.add_space_dims(2);
  ph2.minimize();

  bool ok = (ph1 == ph2);

  print_cons(ph1, "*** ph1 after add_space_dims ***");
  print_gens(ph2, "*** ph2 after add_space_dims ***");

  return ok;
}

bool
test09() {
  Var x(0);

  Poly ph1(1, Topol::NNC);

  ph1.add_con(x > -3);
  ph1.add_con(x < 3);

  print_cons(ph1, "*** ph1 ***");

  ph1.add_space_dims(2, true);

  Gens gs;
  gs.push_back(point());
  gs.push_back(closure_point(-3*x));
  gs.push_back(closure_point(3*x));

  Poly ph2(1, Topol::NNC, Spec_Elem::EMPTY);
  ph2.add_gens(gs);

  print_gens(ph2, "*** ph2 ***");

  ph2.add_space_dims(2, true);

  bool ok = (ph1 == ph2);

  print_cons(ph1, "*** ph1 after add_space_dims ***");
  print_gens(ph2, "*** ph2 after add_space_dims ***");

  return ok;
}

bool
test10() {
  Var A(0);

  Poly ph(1, Topol::NNC, Spec_Elem::EMPTY);
  ph.add_gen(point(A));
  ph.add_gen(closure_point());
  ph.add_gen(closure_point(3*A));

  print_gens(ph, "*** ph ***");

  ph.add_space_dims(1);

  Poly known_result(2, Topol::NNC);
  known_result.add_con(A > 0);
  known_result.add_con(A < 3);

  bool ok = (ph == known_result);

  print_cons(ph, "*** after ph.add_space_dims(1) ***");
  print_gens(ph, "*** after ph.add_space_dims(1) ***");

  return ok;
}

bool
test11() {
  Var A(0);
  Var B(1);

  Gens gs1;
  gs1.push_back(point(A));
  Poly ph(1, Topol::NNC, Spec_Elem::EMPTY);
  ph.add_gens(gs1);

  Gens gs2;
  gs2.push_back(point(A));
  gs2.push_back(ray(A));
  ph.add_gens(gs2);

  print_gens(ph, "*** ph ***");

  ph.add_space_dims(1);

  Gens known_gs;
  known_gs.push_back(point(A));
  known_gs.push_back(line(B));
  known_gs.push_back(ray(A));
  Poly known_result(2, Topol::NNC, Spec_Elem::EMPTY);
  known_result.add_gens(known_gs);

  bool ok = (ph == known_result);

  print_gens(ph, "*** ph ***");

  return ok;
}

bool
test12() {
  Var A(0);
  Var B(1);

  Gens gs1;
  gs1.push_back(point(A + B));
  gs1.push_back(closure_point());
  gs1.push_back(ray(A));
  gs1.push_back(ray(B));
  Poly ph1(2, Topol::NNC, Spec_Elem::EMPTY);
  ph1.add_gens(gs1);

  ph1.minimize();

  print_cons(ph1, "*** ph1 ***");
  print_gens(ph1, "*** ph1 ***");

  ph1.add_space_dims(1);

  Poly known_result(3, Topol::NNC);
  known_result.add_con(A > 0);
  known_result.add_con(B > 0);

  bool ok = (ph1 == known_result);

  print_cons(ph1, "*** after ph1.add_space_dims(1) ***");
  print_gens(ph1, "*** after ph1.add_space_dims(1) ***");

  return ok;
}

bool
test13() {
  Var A(0);
  Poly ph1(2);
  ph1.add_con(A >= 0);
  ph1.add_con(A <= 2);
  ph1.add_space_dims(1);
  Poly& ph1_alias = ph1;
  ph1_alias = ph1; // self-assignment
  return ph1.check_inv();
}

bool
test14() {
  Var A(0);
  Var B(1);
  Var C(2);

  Gens gs1;
  gs1.push_back(point(A));
  gs1.push_back(point(B));
  gs1.push_back(closure_point());
  gs1.push_back(ray(A));
  gs1.push_back(ray(B));
  Poly ph1(2, Topol::NNC, Spec_Elem::EMPTY);
  ph1.add_gens(gs1);

  ph1.add_space_dims(1, true);

  ph1.add_gen(ray(C));
  ph1.add_gen(point(C));

  Poly known_result(3, Topol::NNC);
  known_result.add_con(A >= 0);
  known_result.add_con(B >= 0);
  known_result.add_con(C >= 0);
  known_result.add_con(A + B + C > 0);

  bool ok = (ph1 == known_result);

  ph1.ascii_dump(nout);

  return ok;
}

bool
test15() {
  Var A(0);
  Var B(1);
  Var C(2);

  Gens gs1;
  gs1.push_back(point(A));
  gs1.push_back(point(B));
  gs1.push_back(closure_point());
  gs1.push_back(ray(A));
  gs1.push_back(ray(B));
  Poly ph1(2, Topol::NNC, Spec_Elem::EMPTY);
  ph1.add_gens(gs1);

  ph1.add_space_dims(1, false);

  ph1.add_con(C >= 0);
  ph1.minimize();

  Poly known_result(3, Topol::NNC);
  known_result.add_con(A >= 0);
  known_result.add_con(B >= 0);
  known_result.add_con(C >= 0);
  known_result.add_con(A + B > 0);

  bool ok = (ph1 == known_result);

  ph1.ascii_dump(nout);

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
END_MAIN
