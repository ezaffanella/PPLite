/* PPLite: a lightweight library for convex polyhedra derived from PPL.
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
  auto sd = 4;

  Poly ph(sd);
  for (auto i = 0; i < sd; ++i) {
    ph.add_con(Var(i) >= 0);
    ph.add_con(Var(i) <= 10 + i);
  }
  ph.add_con(Var(0) + Var(1) >= 1);
  ph.minimize();

  F_Poly fp(sd);
  for (auto i = 0; i < sd; ++i) {
    fp.add_con(Var(i) >= 0);
    fp.add_con(Var(i) <= 10 + i);
  }
  fp.add_con(Var(0) + Var(1) >= 1);

  Poly ph_from_fp = fp.to_poly();

  bool ok = ph_from_fp.check_inv()
    && (ph == ph_from_fp);

  return ok;
}

bool
test02() {
  auto sd = 17;
  auto var = [](char c) { return Var(c - 'A'); };

  // The join in this test comes from a PHAVerLite benchmark.

  Cons cs1 = {
    var('A') >= 0, -var('A') >= -1,
    var('B') >= 1, -var('B') >= -4,
    var('C') >= 1, -var('C') >= -8,
    var('D') >= 1, -var('D') >= -12,
    var('E') == 1,
    var('F') >= 1, -var('F') >= -2,
    var('G') >= 1, -var('G') >= -3,
    var('H') == 0,
    var('I') >= 0, -var('I') >= -4,
    var('J') >= 0, -var('J') >= -8,
    var('K') >= 0, -var('K') >= -12,
    var('L') >= 0,
    4*var('M') >= 1, -var('M') >= -1,
    var('N') == 0,
    var('O') >= 0, -var('O') >= -1,
    var('P') >= 0, -var('P') >= -1,
    var('Q') >= 0, -var('Q') >= -1,
    -var('I') - var('L') + var('M') >= -8
  };

  Cons cs2 = {
    var('A') >= 0, -var('A') >= -1,
    var('B') >= 1, -var('B') >= -4,
    var('C') == 1,
    var('D') >= 1, -var('D') >= -12,
    var('E') == 1,
    var('F') >= 1, -var('F') >= -2,
    var('G') == 1,
    var('H') == 0,
    var('I') >= 0, -var('I') >= -4,
    var('J') >= 0, -var('J') >= -8,
    16*var('K') >= 15, -var('K') >= -12,
    var('L') >= 0, -var('L') >= -6,
    4*var('M') >= 1, -var('M') >= -1,
    var('N') == 0,
    var('O') >= 0, -var('O') >= -1,
    var('P') == 1,
    var('Q') >= 0, -var('Q') >= -1,
    -var('J') - var('L') + var('M') >= -21
  };

  // Compute join using Poly.
  nout << "Time for (minimized) join in Poly: " << std::flush;
  Clock clock;
  Poly ph1(sd);
  ph1.add_cons(cs1);
  Poly ph2(sd);
  ph2.add_cons(cs2);
  ph1.poly_hull_assign(ph2);
  ph1.minimize();
  clock.print_elapsed(nout);
  nout << "\n";

  // Compute join using F_Poly and convert back to Poly.
  nout << "Time for (minimized) join in F_Poly: " << std::flush;
  clock.restart();
  F_Poly fp1(sd);
  fp1.add_cons(cs1);
  F_Poly fp2(sd);
  fp2.add_cons(cs2);
  fp1.poly_hull_assign(fp2);
  fp1.minimize();
  clock.print_elapsed(nout);
  nout << "\n";
  nout << "Time for (minimized) conversion to Poly: " << std::flush;
  clock.restart();
  Poly ph_from_fp1 = fp1.to_poly();
  clock.print_elapsed(nout);
  nout << "\n";

  bool ok = ph_from_fp1.check_inv() && (ph1 == ph_from_fp1);

  print_cons(ph1, "=== ph1 ===");
  print_cons(ph_from_fp1, "=== ph from F_Poly ===");

  return ok;
}

bool test03() {

  F_Poly fph(0, Spec_Elem::UNIVERSE);

  Poly ph = fph.to_poly();

  return ph.check_inv();
}

bool test04() {

  Var A(0); Var B(1);
  F_Poly fph(2, Spec_Elem::UNIVERSE);
  fph.add_con(A + B >= 0);
  fph.minimize();
  fph.add_con(A >= 0);
  fph.add_con(A <= 1);
  fph.add_con(B >= 0);
  fph.add_con(B <= 1);

  (void) fph.normalized_cons();

  print_cons(fph, "=== fph ===");

  bool ok = (fph.num_min_cons() == 4)
    // Note: all itvs dims.
    && (num_rows(fph.impl().factors) == 0);
  return ok;
}

BEGIN_MAIN
  DO_TEST(test01);
#ifdef NDEBUG
  DO_TEST(test02);
#endif
  DO_TEST(test03);
  DO_TEST(test04);
END_MAIN
