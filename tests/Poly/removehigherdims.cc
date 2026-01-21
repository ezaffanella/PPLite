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
  ph1.add_con(A >= 0);
  ph1.add_con(A >= 1);
  ph1.add_con(B >= 0);

  ph1.add_gen(point(2*A));

  print_gens(ph1, "*** ph1 : pending ***");
  print_cons(ph1, "*** ph1 ***");
  print_gens(ph1, "*** ph1 ***");

  Poly ph2(2);
  ph2.add_cons(ph1.copy_cons());

  Poly known_result(2, Spec_Elem::EMPTY);
  known_result.add_gen(point(A));
  known_result.add_gen(ray(A));
  known_result.add_gen(ray(B));

  bool ok = (ph2 == known_result);

  return ok;
}

bool
test02() {
  Var A(0);
  Var B(1);
  Var C(2);

  Poly ph1(3, Topol::NNC);
  ph1.add_con(A >= 1);
  ph1.add_con(B >= 1);
  ph1.add_con(C >= 0);
  ph1.add_con(C <= 4);

  print_cons(ph1, "*** ph1 : pending ***");

  ph1.remove_higher_space_dims(2);

  print_gens(ph1, "*** ph1 : gs_only ***");
  print_cons(ph1, "*** ph1 ***");

  Poly known_result(2, Spec_Elem::EMPTY, Topol::NNC);
  known_result.add_gen(point(A + B));
  known_result.add_gen(ray(A));
  known_result.add_gen(ray(B));

  bool ok = (ph1 == known_result);

  return ok;
}

bool
test03() {
  Var A(0);
  Var B(1);
  Var C(2);

  Poly ph1(3, Topol::NNC);
  ph1.add_con(A >= 0);
  ph1.add_con(B >= 0);
  ph1.add_con(A + B <= 2);

  print_cons(ph1, "*** ph1 : pending ***");
  ph1.remove_higher_space_dims(2);
  print_gens(ph1, "*** ph1 : pending after remove ***");

  ph1.add_gen(closure_point(2*A + 2*B));
  print_gens(ph1, "*** ph1 : pending after add_gen ***");

  print_cons(ph1, "*** ph1 ***");

  Poly known_result(2, Spec_Elem::EMPTY, Topol::NNC);
  known_result.add_gen(point());
  known_result.add_gen(point(2*A));
  known_result.add_gen(point(2*B));
  known_result.add_gen(closure_point(2*A + 2*B));

  bool ok = (ph1 == known_result);

  return ok;
}

bool
test04() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Topol::NNC);
  ph1.add_con(A > 0);
  ph1.add_con(A < 2);
  ph1.add_con(B > 0);
  ph1.add_con(B < 2);
  ph1.add_gen(point(B));

  Poly known_result(2, Spec_Elem::EMPTY, Topol::NNC);
  known_result.add_gen(point(A + B));
  known_result.add_gen(closure_point());
  known_result.add_gen(closure_point(2*A));
  known_result.add_gen(closure_point(2*B));
  known_result.add_gen(closure_point(2*A + 2*B));

  bool ok = (ph1 != known_result);

  return ok;
}

bool
test05() {
  Var A(0);
  Var B(1);
  Var C(2);
  Var D(3);
  Var E(4);

  Poly ph1(5, Topol::NNC);
  ph1.add_con(A > 0);
  ph1.add_con(A <= 1);
  ph1.add_con(B >= 0);
  ph1.add_con(B < 1);
  ph1.add_con(C >= 0);
  ph1.add_con(C < 1);
  ph1.add_con(D >= 0);

  ph1.remove_higher_space_dims(1);

  Poly known_result(1, Spec_Elem::EMPTY, Topol::NNC);
  known_result.add_gen(point(A));
  known_result.add_gen(closure_point());

  bool ok = (ph1 == known_result);

  return ok;
}

bool
test06() {
  Var A(0);
  Var B(1);
  Var C(2);
  Var D(3);
  Var E(4);

  Poly ph1(5, Spec_Elem::EMPTY, Topol::NNC);
  ph1.add_gen(point(A + B + C));
  ph1.add_gen(point(A + C + D + E));
  ph1.add_gen(ray(A + D));
  ph1.add_gen(ray(C + E));
  ph1.add_gen(line(B));

  ph1.remove_higher_space_dims(2);

  Poly known_result(2, Spec_Elem::EMPTY, Topol::NNC);
  known_result.add_gen(point(A + B));
  known_result.add_gen(ray(A));
  known_result.add_gen(line(B));

  bool ok = (ph1 == known_result);
  return ok;
}

bool
test07() {
  Var A(0);
  Var B(1);
  Var C(2);
  Var D(3);
  Var E(4);

  Poly ph1(5, Topol::NNC);
  ph1.add_con(A >= 0);
  ph1.add_con(A <= 1);
  ph1.add_con(B > 0);
  ph1.add_con(B < 1);
  ph1.add_con(C > 0);
  ph1.add_con(C < 1);
  ph1.add_con(D > 0);

  ph1.minimize();
  nout << "*** ph1 before: ***" << endl;
  ph1.ascii_dump(nout);
  ph1.remove_higher_space_dims(1);
  nout << "\n*** ph1 after***" << endl;
  ph1.ascii_dump(nout);

  Poly known_result(1, Spec_Elem::EMPTY, Topol::NNC);
  known_result.add_gen(point());
  known_result.add_gen(point(A));

  ph1.minimize();
  known_result.minimize();
  nout << "\n\n*** ph1 ***" << endl;
  ph1.ascii_dump(nout);
  nout << "\n*** known_result ***" << endl;
  known_result.ascii_dump(nout);

  bool ok = (ph1 == known_result);
  return ok;
}

bool
test08() {
  Var A(0);
  Var B(1);
  Var C(2);
  Var D(3);
  Var E(4);

  Poly ph1(5, Topol::NNC);
  ph1.add_con(A > 0);
  ph1.add_con(A < 2);
  ph1.add_con(B >= 0);
  ph1.add_con(B < 1);
  ph1.add_con(C >= 0);
  ph1.add_con(C < 1);
  ph1.add_con(D >= 0);

  ph1.minimize();
  nout << "*** ph1 before: ***" << endl;
  ph1.ascii_dump(nout);
  ph1.remove_higher_space_dims(1);
  nout << "\n*** ph1 after***" << endl;
  ph1.ascii_dump(nout);

  Poly known_result(1, Spec_Elem::EMPTY, Topol::NNC);
  known_result.add_gen(point(A));
  known_result.add_gen(closure_point());
  known_result.add_gen(closure_point(2*A));

  bool ok = (ph1 == known_result);
  return ok;
}

bool
test09() {
  Var x(0);
  Var y(1);
  Var z(2);
  Var w(3);

  Gens gs;
  gs.push_back(point(x + y + 2*z - w));

  Poly ph(4, Spec_Elem::EMPTY);
  ph.add_gens(gs);

  print_gens(ph, "*** ph ***");

  ph.remove_higher_space_dims(2);

  Gens gs_known_result;
  gs_known_result.push_back(point(x + y));
  Poly known_result(2, Spec_Elem::EMPTY);
  known_result.add_gens(gs_known_result);

  bool ok = (ph == known_result);

  print_gens(ph, "*** after remove_higher_space_dims(2) ***");

  return ok;
}

bool
test10() {
  Poly ph(4, Spec_Elem::EMPTY);

  print_cons(ph, "*** ph ***");

  ph.remove_higher_space_dims(0);

  Poly known_result(0, Spec_Elem::EMPTY);

  bool ok = (ph == known_result);

  print_cons(ph, "*** ph after remove_higher_space_dims(0) ***");

  return ok;
}

bool
test11() {
  Var A(0);

  Poly ph(2);
  ph.add_con(A >= 3);

  print_cons(ph, "*** ph ***");

  Poly known_result(ph);

  ph.remove_higher_space_dims(2);

  bool ok = (ph == known_result);

  print_cons(ph, "*** after ph.remove_higher_space_dims(2) ***");

  return ok;
}

bool
test12() {
  Var A(0);
  Var B(1);

  Poly ph(2, Spec_Elem::EMPTY);
  ph.add_gen(point(2*A + B, 4));

  print_cons(ph, "*** ph ***");

  ph.remove_higher_space_dims(1);

  Poly known_result(1, Spec_Elem::EMPTY);
  known_result.add_gen(point(A, 2));

  bool ok = (ph == known_result);

  print_cons(ph, "*** after ph.remove_higher_space_dims(1) ***");

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
END_MAIN
