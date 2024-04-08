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

// Direct encoding for strict sk ineqs.
bool
test01() {
  Var x(0);
  Var y(1);

  Poly ph1(2, Topol::NNC);
  ph1.add_con(x >= 0);
  ph1.add_con(y >= 0);
  ph1.add_con(x + y <= 4);

  Poly ph2(2, Topol::NNC);
  ph2.add_con(x >= 1);
  ph2.add_con(y >= 1);
  ph2.add_con(x + y < 4);

  print_cons(ph1, "=== ph1 ===");
  print_cons(ph2, "=== ph2 ===");

  ph1.widening_assign(ph2);

  print_cons(ph1, "=== ph1.widening_assign(ph2) ===");

  Poly kr(2, Topol::NNC);
  // PPL widening is less precise.
  kr.add_con(x + y <= 4);

  return (ph1 == kr);
}

// Combinatorial encoding for ns strict ineqs
// goes beyond their different geometric repr.
// ns becomes ns as its support is a subset of stable cons.
bool
test02() {
  Var x(0);
  Var y(1);

  Poly ph1(2, Topol::NNC);
  ph1.add_con(x >= 0);
  ph1.add_con(y >= 0);
  ph1.add_con(x + 2*y > 0);

  Poly ph2(2, Topol::NNC);
  ph2.add_con(x >= 0);
  ph2.add_con(y >= 0);
  ph2.add_con(x + y > 0);

  print_cons(ph1, "=== ph1 ===");
  print_cons(ph2, "=== ph2 ===");

  ph1.widening_assign(ph2);

  print_cons(ph1, "=== ph1.widening_assign(ph2) ===");

  Poly kr(2, Topol::NNC);
  kr.add_con(x >= 0);
  kr.add_con(y >= 0);
  // PPL widening is less precise.
  kr.add_con(x + y > 0);

  return (ph1 == kr);
}

// This is a variant of test02, using a single space dimension.
// The name matches the corresponding test written for PPL.
bool
test02_ter() {
  Var x(0);

  Poly ph1(1, Topol::NNC);
  ph1.add_con(x >= 0);
  ph1.add_con(x < 5);

  Poly ph2(1, Topol::NNC);
  ph2.add_con(x >= 0);
  ph2.add_con(x < 1);
  ph2.add_gen(closure_point(5*x));

  print_cons(ph1, "=== ph1 ===");
  print_cons(ph2, "=== ph2 ===");

  ph1.widening_assign(ph2);

  print_cons(ph1, "=== ph1.widening_assign(ph2) ===");

  Poly kr(1, Topol::NNC);
  kr.add_con(x >= 0);
  // PPL widening is less precise.
  kr.add_con(x < 5);

  return (ph1 == kr);
}

// A sk constraint is stable only if saturates a face
// cut by a *sk* constraint.
// No need for the extended version.
bool
test03() {
  Var x(0);
  Var y(1);

  Poly ph1(2, Topol::NNC);
  ph1.add_con(x >= 0);
  ph1.add_con(y >= 0);
  ph1.add_con(x + y >= 4);
  Poly ph2(2, Topol::NNC);
  ph2.add_con(x >= 2);
  ph2.add_con(y >= 2);
  ph2.add_con(x + y > 4);

  print_cons(ph1, "=== ph1 ===");
  print_cons(ph2, "=== ph2 ===");

  ph1.widening_assign(ph2);

  print_cons(ph1, "=== ph1.widening_assign(ph2) ===");

  Poly kr(2, Topol::NNC);
  // As precise as PPL.
  // kr.add_con(x + y >= 4);

  return (ph1 == kr);
}

// A ns constraint is not stable even if saturates a face
// cut by a *ns* constraint (with a different support)
// because its support is not recognized as stable.
// Without the extended version,
// the result is *less* precise than PPL's.
bool
test04() {
  Var x(0);
  Var y(1);

  Poly ph1(2, Topol::NNC);
  ph1.add_con(x >= 0);
  ph1.add_con(y >= 0);
  ph1.add_con(x + y > 0);
  ph1.minimize();

  Poly ph2(2, Topol::NNC);
  ph2.add_con(3*x - y >= 0);
  ph2.add_con(x - 3*y <= 0);
  ph2.add_con(x + 2*y > 0);
  ph2.minimize();

  if (not (ph1.check_inv() && ph2.check_inv() && ph1.contains(ph2)))
    return false;

  print_cons(ph1, "=== ph1 ===");
  print_cons(ph2, "=== ph2 ===");

  ph1.widening_assign(ph2);

  print_cons(ph1, "=== ph1.widening_assign(ph2) ===");

  if (not (ph1.check_inv() && ph1.contains(ph1) && ph1.contains(ph2)))
    return false;

  Poly knres(2, Topol::NNC);
  // Less precise than PPL.
  // knres.add_con(x + y > 0);

  return (ph1 == knres);
}

// Test for efc
bool
test05() {
  Var x(0);
  Var y(1);

  Poly ph1(2, Topol::NNC);
  ph1.add_con(x >= 0);
  ph1.add_con(y >= 0);
  ph1.add_con(x <= 4);
  ph1.add_con(y <= 4);
  ph1.minimize();

  Poly ph2(2, Topol::NNC);
  ph2.add_con(x > 0);
  ph2.add_con(y >= 0);
  ph2.add_con(x <= 4);
  ph2.add_con(y <= 4);
  ph2.minimize();

  print_cons(ph1, "=== ph1 ===");
  print_cons(ph2, "=== ph2 ===");

  ph1.widening_assign(ph2);

  print_cons(ph1, "=== ph1.widening_assign(ph2) ===");

  if (not (ph1.check_inv() && ph1.contains(ph1) && ph1.contains(ph2)))
    return false;

  Poly knres(2, Topol::NNC);
  knres.add_con(x >= 0);
  knres.add_con(y >= 0);
  knres.add_con(x <= 4);
  knres.add_con(y <= 4);

  return (ph1 == knres);
}

// Add a ns if it is a subset of stable constraint,
// even if there is no corresponding constraint.
// ph2 \stdwiden ph1 = ph2
bool
test06() {
  Var x(0);
  Var y(1);
  Var z(2);

  Poly ph1(3, Topol::NNC, Spec_Elem::EMPTY);
  ph1.add_gen(point(x));
  ph1.add_gen(point(y));
  ph1.add_gen(point(z));
  ph1.add_gen(closure_point());
  ph1.minimize();

  Poly ph2(3, Topol::NNC, Spec_Elem::EMPTY);
  ph2.add_gen(point(x));
  ph2.add_gen(point(y));
  ph2.add_gen(closure_point(z));
  ph2.add_gen(closure_point());
  ph2.minimize();

  print_cons(ph1, "=== ph1 ===");
  print_cons(ph2, "=== ph2 ===");

  Poly knres(ph1);

  ph1.widening_assign(ph2);

  print_cons(ph1, "=== ph1.widening_assign(ph2) ===");

  if (not (ph1.check_inv() && ph1.contains(ph1) && ph1.contains(ph2)))
    return false;

  return (ph1 == knres);
}

// Sas example.
bool
test07() {
  Var x(0);
  Var y(1);

  Poly ph1(2, Topol::NNC);
  ph1.add_con(x >= 0);
  ph1.add_con(x <= 2);
  ph1.add_con(y >= 0);
  ph1.add_con(y <= 2);
  ph1.add_con(x + y > 0);
  ph1.minimize();

  Poly ph2(ph1);
  ph2.add_con(x < 2);
  ph2.add_con(x + 2*y <= 4);
  ph2.minimize();

  print_cons(ph1, "=== ph1 ===");
  print_cons(ph2, "=== ph2 ===");

  Poly knres(2, Topol::NNC);
  knres.add_con(x <= 2);
  knres.add_con(x >= 0);
  knres.add_con(y >= 0);
  knres.add_con(x + y > 0);

  ph1.widening_assign(ph2);

  print_cons(ph1, "=== ph1.widening_assign(ph2) ===");

  if (not (ph1.check_inv() && ph1.contains(ph1) && ph1.contains(ph2)))
    return false;

  return (ph1 == knres);
}

bool
test07_bis() {
  Var x(0);
  Var y(1);

  Poly ph1(2, Topol::NNC);
  ph1.add_con(x >= 0);
  ph1.add_con(x <= 2);
  ph1.add_con(y >= 0);
  ph1.add_con(y <= 2);
  ph1.add_con(x + y > 0);
  ph1.minimize();

  Poly ph2(ph1);
  ph2.add_con(x > 0);
  ph2.add_con(x + y <= 2);
  ph2.minimize();

  print_cons(ph1, "=== ph1 ===");
  print_cons(ph2, "=== ph2 ===");

  Poly knres(2, Topol::NNC);
  knres.add_con(x >= 0);
  knres.add_con(y >= 0);
  knres.add_con(x + y > 0);

  ph1.widening_assign(ph2);

  print_cons(ph1, "=== ph1.widening_assign(ph2) ===");

  if (not (ph1.check_inv() && ph1.contains(ph1) && ph1.contains(ph2)))
    return false;

  return (ph1 == knres);
}

bool
test08() {
  Var x(0);
  Var y(1);

  Poly ph1(2, Topol::NNC);
  ph1.add_con(x >= 0);
  ph1.add_con(x <= 2);
  ph1.add_con(y >= 0);
  ph1.add_con(y <= 2);
  ph1.add_con(x + y > 0);
  ph1.minimize();

  Poly ph2(ph1);
  ph2.add_con(x - y > 0);
  ph2.minimize();

  print_cons(ph1, "=== ph1 ===");
  print_cons(ph2, "=== ph2 ===");

  Poly knres(2, Topol::NNC);
  knres.add_con(x <= 2);
  knres.add_con(y >= 0);

  ph1.widening_assign(ph2);

  print_cons(ph1, "=== ph1.widening_assign(ph2) ===");

  if (not (ph1.check_inv() && ph1.contains(ph1) && ph1.contains(ph2)))
    return false;

  return (ph1 == knres);
}

bool
test09() {
  Var x0(0);
  Var x1(1);
  Var x2(2);
  Var x3(3);
  Var x4(4);

  // Mapping variable from gcc to clang order.
  Dims map(5);
  map[0] = 0;
  map[1] = 1;
  map[2] = 3;
  map[3] = 4;
  map[4] = 2;

  Poly ph1(5, Topol::NNC);
  ph1.add_con(-x3 > 0);
  ph1.add_con(-x1 >= -1);
  ph1.add_con(x1 >= 0);
  ph1.add_con(-x0 >= -19);
  Poly ph2(ph1);
  ph2.map_space_dims(map);

  Poly ph1_check(ph1);
  ph1_check.map_space_dims(map);
  bool ok = (ph1_check == ph2);
  if (!ok)
    return false;

  print_cons(ph1, "=== ph1 ===");
  print_cons(ph2, "=== ph2 ===");

  Poly other1(5, Topol::NNC);
  other1.add_con(x0 == 1);
  other1.add_con(x1 == 0);
  other1.add_con(x3 - 10*x4 == 0);
  other1.add_con(-x3 > 0);
  Poly other2(5, Topol::NNC);
  other2.add_con(x0 == 1);
  other2.add_con(x1  == 0);
  other2.add_con(10*x2 - x4 == 0);
  other2.add_con(-x2 > 0);
  // From the simple mapping, it should be:
  // other2.add_con(-x4 > 0);

  // Stressing that they are different representations
  // for the *same* polyhedron, with mapped variables:
  // *syntactic* differences are only due to equalities.
  Poly other1_check(other1);
  other1_check.map_space_dims(map);
  ok = (other1_check == other2);
  if (!ok)
    return false;

  print_cons(other1, "=== other1 ===");
  print_cons(other2, "=== other2 ===");

  ph1.widening_assign(other1);
  print_cons(ph1, "=== ph1.widen(other1) ===");

  ph2.widening_assign(other2);
  print_cons(ph2, "=== ph2.widen(other2) ===");

  ph1.map_space_dims(map);
  ok = (ph1 == ph2);

  return ok;
}

bool
test10() {
  Var x0(0);
  Var x1(1);
  Var x2(2);
  Var x3(3);
  Var x4(4);

  Poly ph1(5, Topol::NNC);
  ph1.add_con(-x3 > 0);
  ph1.add_con(-x1 >= -1);
  ph1.add_con(x1 >= 0);
  ph1.add_con(-x0 >= -19);
  Poly ph2(ph1);

  print_cons(ph1, "=== ph1 ===");
  print_cons(ph2, "=== ph2 ===");

  Poly other1(5, Topol::NNC);
  other1.add_con(x0 == 1);
  other1.add_con(x1 == 0);
  other1.add_con(x3 - 10*x4 == 0);
  Poly other2(other1);

  other1.add_con(-x3 > 0);
  other2.add_con(-x4 > 0);

  // Stressing that they are different representations
  // for the *same* polyhedron,
  // *syntactic* differences only due to equalities.
  bool ok = (other1 == other2);
  if (!ok)
    return false;

  print_cons(other1, "=== other1 ===");
  print_cons(other2, "=== other2 ===");

  // Here ph1 == ph2 and other1 == other2.

  ph1.widening_assign(other1);
  print_cons(ph1, "=== ph1.widen(other1) ===");

  ph2.widening_assign(other2);
  print_cons(ph2, "=== ph2.widen(other2) ===");

  ok = (ph1 == ph2);

  return ok;
}

// Simplified test10 with 2 dimensions.
bool
test11() {
  Var x0(0);
  Var x1(1);

  Poly ph1(2, Topol::NNC);
  ph1.add_con(x0 < 0);
  Poly ph2(ph1);

  print_cons(ph1, "=== ph1 ===");
  print_cons(ph2, "=== ph2 ===");

  Poly other1(2, Topol::NNC);
  other1.add_con(x0 - 2*x1 == 0);
  Poly other2(other1);

  other1.add_con(x0 < 0);
  other2.add_con(x1 < 0);

  // Stressing that they are different representations
  // for the *same* polyhedron,
  // *syntactic* differences only due to equalities.
  bool ok = (other1 == other2);
  if (!ok)
    return false;

  print_cons(other1, "=== other1 ===");
  print_cons(other2, "=== other2 ===");

  // Here ph1 == ph2 and other1 == other2.

  ph1.widening_assign(other1);
  print_cons(ph1, "=== ph1.widen(other1) ===");

  ph2.widening_assign(other2);
  print_cons(ph2, "=== ph2.widen(other2) ===");

  ok = (ph1 == ph2);

  return ok;
}

// Example 6 (Figure 4) in sas18.
bool
test12() {
  Var x(0);
  Var y(1);

  Poly ph1(2, Topol::NNC);
  ph1.add_con(0 <= x);
  ph1.add_con(0 <= y);
  ph1.add_con(y < 1);
  ph1.add_con(x + y <= 1);

  Poly ph2(2, Topol::NNC);
  ph2.add_con(0 <= x);
  ph2.add_con(0 <= y);
  ph2.add_con(y < 1);
  ph2.add_con(x + 2*y <= 2);

  ph1.minimize();
  ph2.minimize();

  print_cons(ph1, "=== ph1 ===");
  print_cons(ph2, "=== ph2 ===");

  // Stressing that the inclusion hypothesis holds.
  bool ok = ph2.contains(ph1);
  if (!ok)
    return false;

  ph2.widening_assign(ph1);
  ph2.minimize();
  print_cons(ph2, "=== ph2.widen(ph1) ===");

  Poly kr(2, Topol::NNC);
  kr.add_con(0 <= x);
  kr.add_con(0 <= y);
  kr.minimize();
  print_cons(kr, "=== known result ===");

  return (ph2 == kr);
}

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
  DO_TEST(test02_ter);
  DO_TEST(test03);
  DO_TEST(test04);
  DO_TEST(test05);
  DO_TEST(test06);
  DO_TEST(test07);
  DO_TEST(test07_bis);
  DO_TEST(test08);
  DO_TEST(test09);
  DO_TEST(test10);
  DO_TEST(test11);
  DO_TEST(test12);
END_MAIN
