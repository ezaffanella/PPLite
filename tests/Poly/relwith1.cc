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

using namespace pplite::IO_Operators;

bool
test01() {
  Var x(0);
  Var y(1);

  Cons cs;
  cs.push_back(2*x - y >= 0);
  cs.push_back(y >= 0);
  Poly ph(2);
  ph.add_cons(cs);
  print_cons(ph, "*** ph ***");

  Gen g = ray(x + y);
  nout << g << endl;

  Poly_Gen_Rel rel = ph.relation_with(g);

  Poly_Gen_Rel known_rel = Poly_Gen_Rel::subsumes();
  bool ok = (rel == known_rel);

  nout << "ph.relation_with(r(A + B)) == " << rel << endl;

  return ok;
}

bool
test02() {
  Var x(0);

  Poly ph(2, Spec_Elem::EMPTY);
  print_cons(ph, "*** ph ***");

  Gen g = point(x);
  nout << g << endl;

  Poly_Gen_Rel rel = ph.relation_with(g);

  Poly_Gen_Rel known_rel = Poly_Gen_Rel::nothing();
  bool ok = (rel == known_rel);

  nout << "ph.relation_with(v(A)) == " << rel << endl;

  return ok;
}

bool
test03() {
  Poly ph;
  print_cons(ph, "*** ph ***");

  Gen g = point();
  nout << g << endl;

  Poly_Gen_Rel rel = ph.relation_with(g);

  Poly_Gen_Rel known_rel = Poly_Gen_Rel::subsumes();
  bool ok = (rel == known_rel);

  nout << "ph.relation_with(v()) == " << rel << endl;

  return ok;
}

bool
test04() {
  Var x(0);
  Var y(1);

  Gens gs;
  gs.push_back(point());
  gs.push_back(ray(y));
  gs.push_back(line(x));
  Poly ph(2, Spec_Elem::EMPTY);
  ph.add_gens(gs);
  print_gens(ph, "*** ph ***");

  Gen g = point(x + y);
  nout << g << endl;

  Poly_Gen_Rel rel = ph.relation_with(g);

  Poly_Gen_Rel known_rel = Poly_Gen_Rel::subsumes();
  bool ok = (rel == known_rel);

  nout << "ph.relation_with(v(A + B)) == " << rel << endl;

  return ok;
}

bool
test05() {
  Var A(0);
  Var B(1);

  Gens gs;
  gs.push_back(point(1*A + 1*B));
  Poly ph(2, Spec_Elem::EMPTY);
  ph.add_gens(gs);
  print_gens(ph, "*** ph ***");

  Poly_Con_Rel rel = ph.relation_with(A >= 0);

  Poly_Con_Rel known_rel = Poly_Con_Rel::is_included();
  bool ok = (rel == known_rel);

  nout << "ph.relation_with(A >= 0) == " << rel << endl;

  return ok;
}

bool
test06() {
  Var y(1);

  Poly ph(2, Spec_Elem::EMPTY);
  print_gens(ph, "*** ph ***");

  Con c(y >= 0);
  nout << c << endl;

  Poly_Con_Rel rel = ph.relation_with(c);

  Poly_Con_Rel known_rel = Poly_Con_Rel::saturates()
    && Poly_Con_Rel::is_included()
    && Poly_Con_Rel::is_disjoint();
  bool ok = (rel == known_rel);

  nout << "ph.relation_with(c) == " << rel << endl;

  return ok;
}

bool
test07() {
  Poly_Con_Rel rel = Poly_Con_Rel::nothing();
  Poly_Con_Rel known_result = Poly_Con_Rel::nothing();

  Poly ph;
  print_gens(ph, "*** ph ***");

  // A false inequality constraint.
  Con c_false1(-1 >= Linear_Expr());
  nout << c_false1 << endl;

  rel = ph.relation_with(c_false1);

  nout << "ph.relation_with(c_false1) == " << rel << endl;

  known_result = Poly_Con_Rel::is_disjoint();
  bool ok = (rel == known_result);

  // A false equality constraint.
  Con c_false2(Linear_Expr() == -5);
  nout << c_false2 << endl;

  rel = ph.relation_with(c_false2);

  nout << "ph.relation_with(c_false2) == " << rel << endl;

  known_result = Poly_Con_Rel::is_disjoint();
  bool ok1 = (rel == known_result);

  // A saturated inequality.
  Con c_saturated1(Linear_Expr() >= 0);
  nout << c_saturated1 << endl;

  rel = ph.relation_with(c_saturated1);

  nout << "ph.relation_with(c_saturated1) == " << rel << endl;

  known_result = Poly_Con_Rel::saturates()
    && Poly_Con_Rel::is_included();
  bool ok2 = (rel == known_result);

  // A saturated equality.
  Con c_saturated2(Linear_Expr() == 0);
  nout << c_saturated2 << endl;

  rel = ph.relation_with(c_saturated2);

  nout << "ph.relation_with(c_saturated2) == " << rel << endl;

  known_result = Poly_Con_Rel::saturates()
    && Poly_Con_Rel::is_included();
  bool ok3 = (rel == known_result);

  // A satisfied inequality which is not saturated.
  Con c_satisfied(Linear_Expr() >= -2);
  nout << c_satisfied << endl;

  rel = ph.relation_with(c_satisfied);

  nout << "ph.relation_with(c_satisfied) == " << rel << endl;

  known_result = Poly_Con_Rel::is_included();
  bool ok4 = (rel == known_result);

  return ok && ok1 && ok2 && ok3 && ok4;
}

bool
test08() {
  Var x(0);
  Var y(1);

  Cons cs;
  cs.push_back(x + y >= 1);
  cs.push_back(y >= 5);

  Poly ph(2);
  ph.add_cons(cs);

  print_gens(ph, "*** ph ***");

  // An equality constraint non-intersecting the polyhedron.
  Con c(y == -1);
  nout << c << endl;

  Poly_Con_Rel rel = ph.relation_with(c);

  nout << "ph.relation_with(c) == " << rel << endl;

  Poly_Con_Rel known_result = Poly_Con_Rel::is_disjoint();
  return rel == known_result;
}

bool
test09() {
  // The zero-dim universe polyhedron.
  Poly ph;
  Poly_Con_Rel rel = ph.relation_with(Linear_Expr() > 0);

  print_gens(ph, "*** ph ***");
  nout << "ph.relation_with(0 > 0) == " << rel << endl;

  Poly_Con_Rel known_result = Poly_Con_Rel::saturates()
    && Poly_Con_Rel::is_disjoint();

  return rel == known_result;
}

bool
test10() {
  // The zero-dim universe polyhedron.
  Poly ph;
  Poly_Con_Rel rel = ph.relation_with(Linear_Expr() > 1);

  print_gens(ph, "*** ph ***");
  nout << "ph.relation_with(0 > 1) == " << rel << endl;

  Poly_Con_Rel known_result = Poly_Con_Rel::is_disjoint();

  return rel == known_result;
}

bool
test11() {
  // The zero-dim universe polyhedron.
  Poly ph;
  Poly_Con_Rel rel = ph.relation_with(Linear_Expr() < 1);

  print_gens(ph, "*** ph ***");
  nout << "ph.relation_with(1 > 0) == " << rel << endl;

  Poly_Con_Rel known_result = Poly_Con_Rel::is_included();

  return rel == known_result;
}

bool
test12() {
  // An empty polyhedron.
  Poly ph(1);
  ph.add_con(Linear_Expr() >= 1);
  Var A(0);
  Poly_Con_Rel rel = ph.relation_with(A > 0);

  print_gens(ph, "*** ph ***");
  nout << "ph.relation_with(A > 0) == " << rel << endl;

  Poly_Con_Rel known_result = Poly_Con_Rel::saturates()
    && Poly_Con_Rel::is_included()
    && Poly_Con_Rel::is_disjoint();

  return rel == known_result;
}

bool
test13() {
  Var A(0);
  Var B(1);
  Cons cs;
  cs.push_back(A + B == 3);
  Poly ph(2);
  ph.add_cons(cs);

  Poly_Con_Rel rel = ph.relation_with(A + B > 3);

  print_gens(ph, "*** ph ***");
  nout << "ph.relation_with(A + B > 3) == " << rel << endl;

  Poly_Con_Rel known_result = Poly_Con_Rel::saturates()
    && Poly_Con_Rel::is_disjoint();

  return rel == known_result;
}

bool
test14() {
  Var A(0);
  Var B(1);
  Cons cs;
  cs.push_back(A + B <= 3);
  Poly ph(2);
  ph.add_cons(cs);

  Poly_Con_Rel rel = ph.relation_with(A + B > 3);

  print_gens(ph, "*** ph ***");
  nout << "ph.relation_with(A + B > 3) == " << rel << endl;

  Poly_Con_Rel known_result = Poly_Con_Rel::is_disjoint();

  return rel == known_result;
}

bool
test15() {
  Var A(0);
  Var B(1);
  Cons cs;
  cs.push_back(A >= 1);
  cs.push_back(B >= 0);
  cs.push_back(A + B <= 3);
  Poly ph(2);
  ph.add_cons(cs);

  Poly_Con_Rel rel = ph.relation_with(A + 2*B < 10);

  print_gens(ph, "*** ph ***");
  nout << "ph.relation_with(A + 2*B < 10) == " << rel << endl;

  Poly_Con_Rel known_result = Poly_Con_Rel::is_included();

  return rel == known_result;
}

bool
test16() {
  Var A(0);
  Var B(1);
  Cons cs;
  cs.push_back(A >= 1);
  cs.push_back(B >= 0);
  cs.push_back(A + B <= 3);
  Poly ph(2);
  ph.add_cons(cs);

  Poly_Con_Rel rel = ph.relation_with(A + B > 1);

  print_gens(ph, "*** ph ***");
  nout << "ph.relation_with(A + B > 1) == " << rel << endl;

  Poly_Con_Rel known_result = Poly_Con_Rel::strictly_intersects();

  return rel == known_result;
}

bool
test17() {
  Var A(0);

  Poly ph(2);
  ph.add_con(A == 0);

  Poly_Gen_Rel rel = ph.relation_with(point(2*A));

  print_gens(ph, "*** ph ***");
  nout << "ph.relation_with(point(2*A)) == " << rel << endl;

  Poly_Gen_Rel known_result = Poly_Gen_Rel::nothing();

  return rel == known_result;
}

bool
test18() {
  Var A(0);
  Var B(1);

  Gens gs;
  gs.push_back(point(A + 0*B));
  gs.push_back(point(3*A));
  Poly ph(2, Spec_Elem::EMPTY);
  ph.add_gens(gs);

  print_gens(ph, "*** ph ***");

  Poly_Con_Rel rel = ph.relation_with(B == 0);

  Poly_Con_Rel known_rel = Poly_Con_Rel::saturates()
    && Poly_Con_Rel::is_included();

  bool ok = (rel == known_rel);

  nout << "ph.relation_with(B == 0) == " << rel << endl;

  return ok;
}

bool
test19() {
  Var A(0);

  Poly ph(2);
  ph.add_con(A >= 0);

  Poly_Gen_Rel rel = ph.relation_with(ray(-A));

  print_cons(ph, "*** ph ***");
  nout << "ph.relation_with(ray(-A)) == " << rel << endl;

  Poly_Gen_Rel known_result = Poly_Gen_Rel::nothing();
  return rel == known_result;
}

bool
test20() {
  Var A(0);

  Poly ph(2);
  ph.add_con(A >= 0);

  Poly_Gen_Rel rel = ph.relation_with(line(A));

  print_cons(ph, "*** ph ***");
  nout << "ph.relation_with(line(A)) == " << rel << endl;

  Poly_Gen_Rel known_result = Poly_Gen_Rel::nothing();
  return rel == known_result;
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
END_MAIN
