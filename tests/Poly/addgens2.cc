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

  Poly ph(2, Topol::NNC);
  Cons cs;
  cs.push_back(A > 0);
  cs.push_back(B >= 0);
  cs.push_back(B <= 1);
  ph.add_cons(cs);
  ph.minimize();

  cs = ph.copy_cons();
  auto strict_before = std::count_if(cs.begin(), cs.end(),
                                     std::mem_fn(&Con::is_strict_inequality));
  print_cons(ph, "=== ph before adding closure point ===");

  ph.add_gen(closure_point(2*B));
  ph.minimize();

  cs = ph.copy_cons();
  auto strict_after = std::count_if(cs.begin(), cs.end(),
                                    std::mem_fn(&Con::is_strict_inequality));
  print_cons(ph, "=== ph after adding closure point ===");

  Poly kr(2, Topol::NNC);
  cs.clear();
  cs.push_back(A > 0);
  cs.push_back(B >= 0);
  cs.push_back(B < 2);
  kr.add_cons(cs);
  kr.minimize();
  print_cons(kr, "=== expected ===");

  return (strict_after == (strict_before + 1)) && (ph == kr);
}

bool
test02() {
  Var A(0);
  Var B(1);

  Poly ph(2, Topol::NNC);
  Cons cs;
  cs.push_back(A > 1);
  cs.push_back(B >= 1);
  cs.push_back(B <= 3);
  ph.add_cons(cs);
  ph.minimize();

  ph.add_gen(closure_point());
  ph.minimize();

  cs = ph.copy_cons();
  int ph_strict = std::count_if(cs.begin(), cs.end(),
                                std::mem_fn(&Con::is_strict_inequality));
  Poly kr(2, Topol::NNC);
  cs.clear();
  cs.push_back(3*A - B > 0);
  cs.push_back(B > 0);
  cs.push_back(B <= 3);
  kr.add_cons(cs);
  kr.minimize();

  cs = kr.copy_cons();
  int kr_strict = std::count_if(cs.begin(), cs.end(),
                                std::mem_fn(&Con::is_strict_inequality));
  return  kr_strict == ph_strict;
}

bool
test04() {
  Var A(0);
  Var B(1);
  Var C(2);

  Poly ph(3, Topol::NNC, Spec_Elem::EMPTY);
  ph.add_gen(point());
  ph.add_gen(point(A));
  ph.add_gen(point(B));
  ph.add_gen(point(A+B));
  ph.add_gen(ray(B));
  ph.minimize();

  ph.add_gen(line(B));
  ph.add_gen(ray(A));
  ph.add_gen(point(-A));

  ph.add_con(A>0);
  ph.add_con(B<=1);
  ph.add_gen(ray(C));
  ph.add_gen(point(2*B));

  ph.add_con(A+B+C <= 5);
  ph.add_gen(point(6*B));
  ph.minimize();

  return ph.check_inv();
}

bool
test05() {
  Var A(0);
  Var B(1);
  Var C(2);

  Poly ph(3, Topol::NNC, Spec_Elem::EMPTY);
  ph.add_gen(point());
  ph.add_gen(point(A));
  ph.add_gen(point(B));
  ph.add_gen(point(A+B));
  ph.add_gen(point(C));
  ph.add_gen(point(C+A));
  ph.add_gen(point(C+B));
  ph.add_gen(point(C+A+B));
  ph.minimize();

  ph.add_gen(ray(A+B));
  ph.add_gen(ray(A));
  ph.add_gen(ray(B));
  ph.add_gen(ray(C));
  ph.minimize();

  ph.add_gen(point(-A));
  ph.add_gen(point(-B));
  ph.add_gen(point(-C));
  ph.add_gen(point(-B-C));
  ph.add_gen(point(-A-B));
  ph.add_gen(point(-A-C));
  ph.add_gen(point(-A-B-C));
  ph.minimize();

  return ph.check_inv();
}

/* Note:
 * The commented stats are to be checked in the opt build
 * since in the debug one many more conversions may be called.
 */
bool
test06() {
  Var A(0);
  Var B(1);
  Var C(2);
  Var D(3);

  Poly ph(4, Topol::NNC, Spec_Elem::EMPTY);
  ph.add_gen(point());
  ph.add_gen(point(A));
  ph.add_gen(point(B));
  ph.minimize();

  // point + efc
  ph.add_gen(point(A+B));
  ph.add_gen(point(2*A + 2*B));
  ph.add_gen(point(2*A));
  ph.add_gen(point(2*B));
  ph.minimize();  // 4/4 - was 4/4

  // ray + efc
  ph.add_gen(ray(B));
  // point + rfc
  ph.add_gen(point(A - B));
  // ph.minimize(); // 2/2 - was 0/2

  // ray + rfc
  ph.add_gen(ray(A));
  // ph.minimize(); // 0/1 - was 0/1

  // ray + positivity
  ph.add_gen(ray(-A + B));
  // ph.minimize(); // 1/1 - was 0/1

  // point + positivity
  ph.add_gen(point(C));
  // ph.minimize(); // 1/1 - was 0/1

  // point + rfc
  ph.add_gen(point(A-B+C));
  // ph.minimize(); // 0/1 - was 0/1

  ph.minimize();  // tot: 4/5 - was 0/5

  // line + rfc
  ph.add_gen(line(C));
  ph.minimize();  // 1/1 - was 0/1

  return ph.check_inv();
}

bool
test07() {
  Cons cs;
  cs.push_back(Var(0) -Var(9) == 0);
  cs.push_back(-Var(2) -7*Var(4) >= -15);
  cs.push_back(Var(0) >= 0);
  cs.push_back(-12*Var(4) -Var(8) >= -24);
  cs.push_back(Var(2) >= 1);
  cs.push_back(Var(1) >= 1);
  cs.push_back(-Var(1) +3*Var(4) >= 2);
  cs.push_back(-4*Var(0) >= -3);
  Poly ph1(13, Topol::NNC);
  ph1.add_cons(cs);

  Gens gs;
  gs.push_back(point(Var(1) + 8*Var(2) + 12*Var(8)));
  gs.push_back(line(Var(3)));
  gs.push_back(line(Var(4)));
  gs.push_back(line(Var(5)));
  gs.push_back(line(Var(6)));
  gs.push_back(line(Var(7)));
  gs.push_back(line(Var(10)));
  gs.push_back(line(Var(11)));
  gs.push_back(line(Var(12)));
  gs.push_back(ray(-Var(8)));
  gs.push_back(point(Var(1) + Var(2) + 12*Var(8)));
  gs.push_back(point(4*Var(1) + Var(2) + 12*Var(8)));
  gs.push_back(point(4*Var(1) + 8*Var(2) + 12*Var(8)));
  gs.push_back(point(Var(1) + 8*Var(2) + 12*Var(8)));
  gs.push_back(point(Var(1) + Var(2) + 12*Var(8)));
  gs.push_back(point(4*Var(1) + Var(2) + 12*Var(8)));
  gs.push_back(point(4*Var(1) + 8*Var(2) + 12*Var(8)));
  ph1.add_gens(gs);

  ph1.minimize();

  return ph1.check_inv();
}

bool
test08() {
  Cons cs;
  cs.push_back(-Var(1) -7*Var(3) >= -15);
  cs.push_back( Var(4) >= 0);
  cs.push_back(-Var(2) -12*Var(3) >= -24);
  cs.push_back(Var(1) >= 1);
  cs.push_back(Var(0) >= 1);
  cs.push_back(-Var(0) +3*Var(3) >= 2);
  cs.push_back(-4*Var(4) >= -3);
  Poly ph1(5, Topol::NNC);
  ph1.add_cons(cs);

  ph1.add_gen(line(Var(3)));

  ph1.minimize();

  return ph1.check_inv();
}

bool
test09() {
  Var A(0);
  Var B(1);
  Var C(2);

  // This example was triggering, when executed in debugging mode,
  // a spurious assertion failure (due to wrong assertion).
  Poly ph(3, Spec_Elem::EMPTY);
  ph.add_gen(point());
  ph.add_gen(point(B));
  ph.add_gen(point(A));
  ph.add_gen(ray(C));
  ph.add_gen(ray(-C));

  // Minimizing here was hiding the assertion failure.
  /* ph.minimize(); */

  ph.add_gen(ray(A));
  ph.minimize();

  Poly kr(3);
  kr.add_cons({ A >= 0, 0 <= B, B <= 1 });

  return (ph == kr);
}

bool
test10() {
  set_default_topology(Topol::NNC);

  Gens gs;
  gs.push_back(point());
  gs.push_back(closure_point(Var(0)));

  Poly ph(1);
  ph.add_gens(gs);
  ph.minimize();

  set_default_topology(Topol::CLOSED);

  return ph.check_inv();
}

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
/* test03 moved into Poly_Impl */
  DO_TEST(test04);
  DO_TEST(test05);
  DO_TEST(test06);
  DO_TEST(test07);
  DO_TEST(test08);
  DO_TEST(test09);
  DO_TEST(test10);
END_MAIN
