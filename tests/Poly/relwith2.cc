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

using namespace pplite::IO_Operators;

bool
test01() {
  Var A(0);
  Var B(1);

  Gens gs1;
  gs1.push_back(point());
  gs1.push_back(line(A + B));
  Poly ph1(2, Spec_Elem::EMPTY);
  ph1.add_gens(gs1);

  Gens gs2;
  gs2.push_back(ray(A));
  gs2.push_back(point(B));
  gs2.push_back(point(-B));
  Poly ph2(2, Spec_Elem::EMPTY);
  ph2.add_gens(gs2);

  Poly_Con_Rel rel1 = ph1.relation_with(A == 0);
  Poly_Con_Rel rel2 = ph2.relation_with(A == 0);

  print_gens(ph1, "*** ph1 ***");
  print_gens(ph2, "*** ph2 ***");
  nout << "ph1.relation_with(A == 0) == " << rel1 << endl;
  nout << "ph2.relation_with(A == 0) == " << rel2 << endl;

  Poly_Con_Rel known_result = Poly_Con_Rel::strictly_intersects();
  return rel1 == known_result && rel2 == known_result;
}

bool
test02() {
  Var A(0);
  Var B(1);

  Gens gs;
  gs.push_back(point(A));
  gs.push_back(line(B));
  Poly ph(2, Spec_Elem::EMPTY);
  ph.add_gens(gs);

  Poly_Con_Rel rel = ph.relation_with(B > 0);

  print_gens(ph, "*** ph ***");
  nout << "ph.relation_with(B > 0) == " << rel << endl;

  Poly_Con_Rel known_result = Poly_Con_Rel::strictly_intersects();
  return rel == known_result;
}

bool
test03() {
  Var A(0);
  Var B(1);

  Poly ph(2);
  ph.add_con(A >= 2);
  ph.add_con(B == 0);

  Poly_Gen_Rel rel = ph.relation_with(ray(A + B));

  Poly_Gen_Rel known_rel = Poly_Gen_Rel::nothing();

  bool ok = (rel == known_rel);

  print_cons(ph, "*** ph ***");
  nout << "ph.relation_with(ray(A + B)) == " << rel << endl;

  return ok;
}

bool
test04() {
  Var A(0);
  Var B(1);

  Poly ph(2, Spec_Elem::EMPTY);
  ph.add_gen(point());
  ph.minimize();
  ph.add_gen(ray(A));
  ph.add_gen(ray(B));

  Poly_Con_Rel rel = ph.relation_with(A == 0);

  Poly_Con_Rel known_rel = Poly_Con_Rel::strictly_intersects();

  bool ok = (rel == known_rel);

  print_cons(ph, "*** ph ***");
  nout << "ph.relation_with(A == 0) == " << rel << endl;

  return ok;
}

bool
test05() {
  Var A(0);
  Var B(1);

  Poly ph(2, Topol::NNC);
  ph.add_con(A - B > 0);
  ph.add_con(B >= 0);

  Poly_Gen_Rel rel1 = ph.relation_with(point(B));
  Poly_Gen_Rel rel2 = ph.relation_with(point(-B));

  print_gens(ph, "*** ph ***");
  nout << "ph.relation_with(point(B)) == " << rel1 << endl;
  nout << "ph.relation_with(point(-B)) == " << rel2 << endl;

  Poly_Gen_Rel known_result = Poly_Gen_Rel::nothing();
  return rel1 == known_result && rel2 == known_result;
}

bool
test06() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::EMPTY, Topol::NNC);
  ph1.add_gen(point(A + B));

  Poly_Con_Rel rel = ph1.relation_with(A - B == 0);
  Poly_Con_Rel known_rel = Poly_Con_Rel::saturates()
    && Poly_Con_Rel::is_included();

  bool ok = (rel == known_rel);

  print_gens(ph1, "*** ph1 ***");
  nout << "ph1.relation_with(A - B == 0) = " << rel << endl;

  return ok;
}

bool
test07() {
  Var A(0);
  Var B(1);

  Gens gs;
  gs.push_back(point());
  gs.push_back(ray(A));
  gs.push_back(line(B));

  Poly ph1(2, Topol::NNC);
  ph1.add_gens(gs);

  Poly_Con_Rel rel = ph1.relation_with(A >= 1);
  Poly_Con_Rel known_rel = Poly_Con_Rel::strictly_intersects();

  bool ok = (rel == known_rel);

  print_gens(ph1, "*** ph1 ***");
  nout << "ph1.relation_with(A >= 1) = " << rel << endl;

  return ok;
}

bool
test08() {
  Var A(0);
  Var B(1);

  Gens gs;
  gs.push_back(line(A));
  gs.push_back(ray(B));
  gs.push_back(point());
  Poly ph1(2, Topol::NNC);
  ph1.add_gens(gs);

  Poly_Con_Rel rel = ph1.relation_with(A > 1);
  Poly_Con_Rel known_rel = Poly_Con_Rel::strictly_intersects();

  bool ok = (rel == known_rel);

  print_gens(ph1, "*** ph1 ***");
  nout << "ph1.relation_with(A > 1) = " << rel << endl;

  return ok;
}

bool
test09() {
  Var A(0);
  Var B(1);

  Gens gs;
  gs.push_back(point(A));
  Poly ph(2, Topol::NNC, Spec_Elem::EMPTY);
  ph.add_gens(gs);

  Poly ph1(2);
  ph1.add_con(A == 1);
  ph1.add_con(B == 1);

  gs = ph.copy_gens();
  Gens gs1;
  std::copy_if(gs.begin(), gs.end(), std::back_inserter(gs1),
               [](const Gen& g) { return !g.is_closure_point(); });
  ph1.add_gens(gs1);

  Poly_Con_Rel rel = ph1.relation_with(B == 1);
  Poly_Con_Rel known_rel = Poly_Con_Rel::strictly_intersects();

  bool ok = (rel == known_rel);

  print_gens(ph1, "*** ph1 ***");
  nout << "ph1.relation_with(B == 1) = " << rel << endl;

  return ok;
}

bool
test10() {
  Var A(0);
  Var B(1);

  Gens gs;
  gs.push_back(ray(B));
  gs.push_back(point(-A));
  Poly ph1(2, Spec_Elem::EMPTY);
  ph1.add_gens(gs);
  ph1.minimize();

  Poly_Con_Rel rel = ph1.relation_with(B <= 0);
  Poly_Con_Rel known_rel = Poly_Con_Rel::strictly_intersects();

  bool ok = (rel == known_rel);

  print_gens(ph1, "*** ph1 ***");
  nout << "ph1.relation_with(B <= 0) = " << rel << endl;

  return ok;
}

bool
test11() {
  Var A(0);
  Var B(1);

  Gens gs;
  gs.push_back(point(A + B));
  gs.push_back(point(-A + B));
  Poly ph1(2, Spec_Elem::EMPTY);
  ph1.add_gens(gs);
  ph1.minimize();

  Poly_Con_Rel rel = ph1.relation_with(A >= 0);
  Poly_Con_Rel known_rel = Poly_Con_Rel::strictly_intersects();

  bool ok = (rel == known_rel);

  print_gens(ph1, "*** ph1 ***");
  nout << "ph1.relation_with(A >= 0) = " << rel << endl;

  return ok;
}

bool
test12() {
  Var A(0);
  Var B(1);

  Gens gs;
  gs.push_back(point(-A));
  gs.push_back(ray(-B));
  gs.push_back(ray(A + B));
  Poly ph(2, Topol::NNC, Spec_Elem::EMPTY);
  ph.add_gens(gs);
  ph.minimize();

  Poly_Con_Rel rel = ph.relation_with(B < 0);
  Poly_Con_Rel known_rel = Poly_Con_Rel::strictly_intersects();

  bool ok = (rel == known_rel);

  print_gens(ph, "*** ph ***");
  nout << "ph.relation_with(B < 0) = " << rel << endl;

  return ok;
}

bool
test13() {
  Var A(0);
  Var B(1);

  Gens gs;
  gs.push_back(ray(A));
  gs.push_back(ray(A + B));
  gs.push_back(point(-B));
  Poly ph1(2, Spec_Elem::EMPTY, Topol::NNC);
  ph1.add_gens(gs);
  ph1.minimize();

  Poly_Con_Rel rel = ph1.relation_with(A < 0);
  Poly_Con_Rel known_rel = Poly_Con_Rel::is_disjoint();

  bool ok = (rel == known_rel);

  print_gens(ph1, "*** ph1 ***");
  nout << "ph1.relation_with(A < 0) = " << rel << endl;

  return ok;
}

bool
test14() {
  Var A(0);
  Var B(1);

  Poly ph(2, Topol::NNC);
  ph.add_con(A == 0);
  ph.add_con(B == 0);

  Poly_Gen_Rel rel = ph.relation_with(closure_point(A));

  print_cons(ph, "*** ph ***");
  nout << "ph.relation_with(line(A)) == " << rel << endl;

  Poly_Gen_Rel known_result = Poly_Gen_Rel::nothing();
  return rel == known_result;
}

bool
test15() {
  Var A(0);
  Var B(1);
  Var C(2);

  Poly ph(3);
  ph.add_con(A >= B);
  ph.add_gen(point());
  ph.minimize();
  ph.add_con(A - B <= -1);

  Poly_Gen_Rel rel = ph.relation_with(line(C));

  print_cons(ph, "*** ph ***");
  nout << "ph.relation_with(line(C)) == " << rel << endl;

  Poly_Gen_Rel known_result = Poly_Gen_Rel::nothing();
  return rel == known_result;
}

bool
test16() {
  Var A(0);

  Poly ph(1);
  ph.add_con(A >= 0);

  Poly_Gen_Rel rel1 = ph.relation_with(ray(A));
  Poly_Gen_Rel rel2 = ph.relation_with(ray(-A));
  return rel1.implies(Poly_Gen_Rel::subsumes())
    && !rel2.implies(Poly_Gen_Rel::subsumes());
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
END_MAIN
