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
  Var A(0), B(1), C(2), D(3);

  Poly ph(4);
  ph.add_con(A >= 0);
  ph.add_con(A <= 10);
  ph.add_con(B >= 0);
  ph.add_con(B <= 10);
  ph.add_con(C >= 0);
  ph.add_con(C <= 10);
  ph.add_con(D >= 0);
  ph.add_con(D <= 10);

  print_cons(ph, "*** ph ***");

  Vars vars = { A, B, C, D };
  Linear_Exprs exprs = { B, C, D, D };
  Integers inhomos = { 0, 1, 2, 3 };
  Integers dens = { 1, 1, 1, 1 };

  ph.parallel_affine_image(vars, exprs, inhomos, dens);

  Poly kr(4);
  kr.add_con(A >= 0);
  kr.add_con(A <= 10);
  kr.add_con(B >= 1);
  kr.add_con(B <= 11);
  kr.add_con(C >= 2);
  kr.add_con(C <= 12);
  kr.add_con(D - C == 1);

  bool ok = (ph == kr);

  print_cons(ph, "*** ph after parallel_affine_image ***");

  return ok;
}

bool
test02() {
  Var A(0), B(1), C(2), D(3);

  Poly ph(4);
  ph.add_con(A >= 0);
  ph.add_con(A <= 10);
  ph.add_con(B >= 0);
  ph.add_con(B <= 10);
  ph.add_con(C >= 0);
  ph.add_con(C <= 10);
  ph.add_con(D >= 0);
  ph.add_con(D <= 10);

  print_cons(ph, "*** ph ***");

  Vars vars = { A, B, C, D };
  Linear_Exprs exprs = { B, C, D, Linear_Expr() };
  Integers inhomos = { 0, 1, 2, 3 };
  Integers dens = { 1, 1, 1, 1 };

  ph.parallel_affine_image(vars, exprs, inhomos, dens);

  Poly kr(4);
  kr.add_con(A >= 0);
  kr.add_con(A <= 10);
  kr.add_con(B >= 1);
  kr.add_con(B <= 11);
  kr.add_con(C >= 2);
  kr.add_con(C <= 12);
  kr.add_con(D == 3);

  bool ok = (ph == kr);

  print_cons(ph, "*** ph after parallel_affine_image ***");

  return ok;
}

bool
test03() {
  Var A(0), B(1), C(2), D(3);

  Poly ph(4);
  ph.add_con(A >= 0);
  ph.add_con(A <= 10);
  ph.add_con(B >= 0);
  ph.add_con(B <= 10);
  ph.add_con(C >= 0);
  ph.add_con(C <= 10);
  ph.add_con(D >= 0);
  ph.add_con(D <= 10);

  print_cons(ph, "*** ph ***");

  Vars vars = { A, C, D, B };
  Linear_Exprs exprs = { B, D, Linear_Expr(), C };
  Integers inhomos = { 0, 2, 3, 1 };
  Integers dens = { 1, 1, 1, 1 };

  ph.parallel_affine_image(vars, exprs, inhomos, dens);

  Poly kr(4);
  kr.add_con(A >= 0);
  kr.add_con(A <= 10);
  kr.add_con(B >= 1);
  kr.add_con(B <= 11);
  kr.add_con(C >= 2);
  kr.add_con(C <= 12);
  kr.add_con(D == 3);

  bool ok = (ph == kr);

  print_cons(ph, "*** ph after parallel_affine_image ***");

  return ok;
}

bool
test04() {
  Var A(0), B(1), C(2), D(3);

  Poly ph(4);
  ph.add_con(A >= 0);
  ph.add_con(A <= 10);
  ph.add_con(B >= 0);
  ph.add_con(B <= 10);
  ph.add_con(C >= 0);
  ph.add_con(C <= 10);
  ph.add_con(D >= 0);
  ph.add_con(D <= 10);

  print_cons(ph, "*** ph ***");

  Vars vars = { A, B, C, D };
  Linear_Exprs exprs = { B, C, D, A };
  Integers inhomos = { 0, 1, 2, 3 };
  Integers dens = { 1, 1, 1, 1 };

  ph.parallel_affine_image(vars, exprs, inhomos, dens);

  Poly kr(4);
  kr.add_con(A >= 0);
  kr.add_con(A <= 10);
  kr.add_con(B >= 1);
  kr.add_con(B <= 11);
  kr.add_con(C >= 2);
  kr.add_con(C <= 12);
  kr.add_con(D >= 3);
  kr.add_con(D <= 13);

  bool ok = (ph == kr);

  print_cons(ph, "*** ph after parallel_affine_image ***");

  return ok;
}

bool
test05() {
  dim_type sd = 4;
  Vars x;
  for (auto i = 0; i != sd; ++i)
    x.push_back(Var(i));

  Poly ph(sd);

  // the discrete transition (in relational style)
  // E - B == 0, F - A == 0, G == 0, H == 0
  Linear_Exprs exprs = { x[1], x[0],  Linear_Expr(), Linear_Expr() };
  Integers inhomos = { 0, 0, 0, 0 };
  Integers dens = { 1, 1, 1, 1 };

  // Implementing the affine image with the ``relational'' approach
  // (i.e., adding sd new dims, adding constraints and then projecting).
  Cons cs_rel;
  for (dim_type i = 0; i != sd; ++i) {
    auto c = (dens[i] * Var(i + sd) - exprs[i] == inhomos[i]);
    cs_rel.push_back(c);
  }

  Poly ph_rel(sd);
  ph_rel.add_space_dims(sd);
  ph_rel.add_cons(cs_rel);
  Dims tbr(sd);
  std::iota(tbr.begin(), tbr.end(), 0);
  ph_rel.remove_space_dims(tbr.begin(), tbr.end());
  ph_rel.minimize();
  print_cons(ph_rel, "*** ph_rel after ***");

  // Reversing the order of parallel assignments should be semantically
  // irrelevant, but it was showing a bug (now corrected).
  std::reverse(x.begin(), x.end());
  std::reverse(exprs.begin(), exprs.end());
  std::reverse(inhomos.begin(), inhomos.end());
  std::reverse(dens.begin(), dens.end());

  Poly ph_par(sd);
  ph_par.parallel_affine_image(x, exprs, inhomos, dens);
  ph_par.minimize();
  print_cons(ph_par, "*** ph_par after ***");

  return (ph_par == ph_rel);
}

bool
test06() {
  dim_type sd = 15;
  Vars x;
  for (auto i = 0; i != sd; ++i)
    x.push_back(Var(i));

  // Init
  // B - G - O = 0,
  // A >= 0, -A >= -5, B >= 0, -B >= -5, C >= 0, -C >= -5, D >= 0, -D >= -5,
  // E >= 0, -E >= -5, F >= 0, -F >= -5, G >= 0, -G >= -5, H >= 0, -H >= -5,
  // I >= 0, -I >= -5, J >= 0, -J >= -5, K >= 0, -K >= -5, L >= 0, -L >= -5,
  // M >= 0, -M >= -5, N >= 0, -N >= -5
  Cons cs;
  for (auto i = 0; i < sd - 1; ++i) {
    cs.push_back(x[i] >= 0);
    cs.push_back(x[i] <= 5);
  }
  cs.push_back(x[1] - x[6] - x[14] == 0);
  Poly ph(sd);
  ph.add_cons(cs);
  // print_cons(ph, "*** ph ***");

  // guard
  // -A + F >= -32, G - H >= -1, -E + G >= -18, K - L >= 3
  Cons guard = {
    - x[0] + x[5] >= -32,
    x[6] - x[7] >= -1,
    - x[4] + x[6] >= -18,
    x[10] - x[11] >= 3
  };
  ph.add_cons(guard);
  ph.minimize();
  // print_cons(ph, "*** ph + guard ***");

  // the discrete transition (in relational style)
  // B - G - D1 = 0, M - C1 = -16, K - B1 = -90, J - A1 = -58, J - Z = -54,
  // I - Y = -35, G - X = -55, N - W = -36, N - V = -38, E - U = -53,
  // D - T = -41, C - S = -41, B - R = -50, F - Q = -72, F - P = -38
  Linear_Exprs exprs = {
    x[5], x[5],  x[1],  x[2],  x[3],
    x[4], x[13], x[13], x[6],  x[8],
    x[9], x[9],  x[10], x[12], x[1] - x[6]
  };
  Integers inhomos = {
    38, 72, 50, 41, 41,
    53, 38, 36, 55, 35,
    54, 58, 90, 16,  0
  };
  Integers dens(sd, 1);

  // Implementing the affine image with the ``relational'' approach
  // (i.e., adding sd new dims, adding constraints and then projecting).
  Cons cs_rel;
  for (dim_type i = 0; i != sd; ++i) {
    auto c = (dens[i] * Var(i + sd) - exprs[i] == inhomos[i]);
    cs_rel.push_back(c);
  }

  // expected result
  // K - L = -4, G - H = 2, C - I - O = -5, A - B = -34, -A >= -43,
  // -C >= -55, -D >= -46, -E >= -46, -F >= -58, -G >= -43, -I >= -60,
  // -J >= -40, -K >= -59, -M >= -95, -N >= -21, N >= 16, M >= 93, K >= 54,
  // J >= 35, I >= 55, G >= 38, F >= 53, E >= 41, D >= 41, C >= 50, A >= 38
  Poly kr(sd);
  Cons kr_cs = {
    x[10] - x[11] == -4,
    x[6] - x[7] == 2,
    x[2] - x[8] - x[14] == -5,
    x[0] - x[1] == -34,
    - x[0] >= -43,
    - x[2] >= -55,
    - x[3] >= -46,
    - x[4] >= -46,
    - x[5] >= -58,
    - x[6] >= -43,
    - x[8] >= -60,
    - x[9] >= -40,
    - x[10] >= -59,
    - x[12] >= -95,
    - x[13] >= -21,
    x[13] >= 16,
    x[12] >= 93,
    x[10] >= 54,
    x[9] >= 35,
    x[8] >= 55,
    x[6] >= 38,
    x[5] >= 53,
    x[4] >= 41,
    x[3] >= 41,
    x[2] >= 50,
    x[0] >= 38
  };
  kr.add_cons(kr_cs);
  kr.minimize();
  // print_cons(kr, "*** expected ***");

  auto ph_rel = ph;
  {
    nout << "Time in relational-style affine image: " << std::flush;
    Clock clock;
    ph_rel.add_space_dims(sd);
    ph_rel.add_cons(cs_rel);
    Dims tbr(sd);
    std::iota(tbr.begin(), tbr.end(), 0);
    ph_rel.remove_space_dims(tbr.begin(), tbr.end());
    ph_rel.minimize();
    clock.print_elapsed(nout);
    nout << "\n";
  }
  // print_cons(ph_rel, "*** ph_rel after ***");

  // Reversing the order of parallel assignments should be semantically
  // irrelevant, but it was showing a bug (now corrected).
  std::reverse(x.begin(), x.end());
  std::reverse(exprs.begin(), exprs.end());
  std::reverse(inhomos.begin(), inhomos.end());
  std::reverse(dens.begin(), dens.end());

  auto ph_par = ph;
  {
    nout << "Time in parallel affine image: " << std::flush;
    Clock clock;
    ph_par.parallel_affine_image(x, exprs, inhomos, dens);
    ph_par.minimize();
    clock.print_elapsed(nout);
    nout << "\n";
  }
  // print_cons(ph_par, "*** ph_par after ***");

  bool ok_par_rel = (ph_par == ph_rel);
  if (!ok_par_rel)
    nout << "\nParallel affine image differs from relational affine image\n";
  bool ok_par = (ph_par == kr);
  if (!ok_par)
    nout << "\nParallel affine image differs from expected\n";
  bool ok_rel = (ph_rel == kr);
  if (!ok_rel)
    nout << "\nRelational affine image differs from expected\n";

  return ok_par_rel && ok_par && ok_rel;
}

bool
test07() {
  dim_type sd = 14;
  Vars x;
  for (auto i = 0; i != sd; ++i)
    x.push_back(Var(i));

  Cons cs;
  for (auto i = 0; i < sd; ++i) {
    cs.push_back(x[i] >= 0);
    cs.push_back(x[i] <= 5);
  }
  Poly ph(sd);
  ph.add_cons(cs);
  ph.minimize();

  Linear_Exprs exprs(sd);
  Integers inhomos(sd);
  std::iota(inhomos.begin(), inhomos.end(), 0);
  Integers dens(sd, 1);

  auto ph1 = ph;
  {
    nout << "Time in parallel affine image: " << std::flush;
    Clock clock;
    ph1.parallel_affine_image(x, exprs, inhomos, dens);
    ph1.minimize();
    clock.print_elapsed(nout);
    nout << "\n";
  }

  auto ph2 = ph;
  {
    nout << "Time in sequence of affine images: " << std::flush;
    Clock clock;
    for (auto i = 0; i < sd; ++i)
      ph2.affine_image(x[i], exprs[i], inhomos[i], dens[i]);
    ph2.minimize();
    clock.print_elapsed(nout);
    nout << "\n";
  }
  return (ph1 == ph2);
}

bool
test08() {
  Var A(0);
  Var B(1);

  Poly ph(2);
  ph.add_con(A >= 0);
  ph.add_con(A <= 10);
  ph.add_con(B >= 0);
  ph.add_con(B <= 10);

  print_cons(ph, "*** ph ***");

  Vars vars = { A, B };
  Linear_Exprs exprs = { A, B };
  Integers inhomos = { 0, 1 };
  Integers dens = { 2, 2 };

  ph.parallel_affine_image(vars, exprs, inhomos, dens);

  Poly kr(2);
  kr.add_con(A >= 0);
  kr.add_con(A <= 5);
  kr.add_con(2*B >= 1);
  kr.add_con(2*B <= 11);

  bool ok = (ph == kr);

  print_cons(ph, "*** ph after parallel_affine_image ***");

  return ok;
}

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
  DO_TEST(test03);
  DO_TEST(test04);
  DO_TEST(test05);
  if (check_exp_eval()) {
    DO_TEST(test06);
  }
#ifdef NDEBUG
  DO_TEST(test07);
#endif
  DO_TEST(test08);
END_MAIN
