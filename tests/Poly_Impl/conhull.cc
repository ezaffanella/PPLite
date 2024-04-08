/* PPLite: a lightweight library for convex polyhedra derived from PPL.
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

namespace {

Poly
create_random_box(dim_type dim) {
  Poly ph(dim);
  for (auto d = 0; d < dim; ++d) {
    ph.add_con((rand() % 20) * Var(d) <= (rand() % 100));
    ph.add_con((rand() % 20) * Var(d) >= -(rand() % 100));
  }
  ph.minimize();
  return ph;
}

} // namespace

bool
test01() {
  dim_type nr = 10;
  dim_type dim = 5;
  if (check_exp_eval()) {
    nr = 30;
    dim = 10;
  }

  std::vector<Poly> phs;

  nout << "\nVector of " << nr << " random "
       << dim << "-dim boxes." << std::endl;
  for (auto i = 0; i < nr; ++i) {
    phs.push_back(create_random_box(dim));
  }

  Poly nary_chull(dim, Spec_Elem::EMPTY);
  {
    Clock clock;
    nary_chull = con_hull(phs.begin(), phs.end());
    nout << "  Time spent in n-ary con hull: ";
    clock.print_elapsed(nout);
    nout << std::endl;
  }

  Poly nary_boxed_chull(dim, Spec_Elem::EMPTY);
  {
    Clock clock;
    nary_boxed_chull = con_hull(phs.begin(), phs.end(), true);
    nout << "  Time spent in n-ary boxed con hull: ";
    clock.print_elapsed(nout);
    nout << std::endl;
  }

  Poly seq_chull(dim, Spec_Elem::EMPTY);
  {
    Clock clock;
    for (auto i = 0; i < nr; ++i) {
      seq_chull.con_hull_assign(phs[i]);
    }
    nout << "  Time spent in sequential con hull: ";
    clock.print_elapsed(nout);
    nout << std::endl;
  }

  bool ok = (nary_boxed_chull == nary_chull
             && nary_chull == seq_chull);
  if (!ok)
    return false;

  ok = std::all_of(phs.begin(), phs.end(), [&nary_chull] (const Poly& ph)
                   { return nary_chull.contains(ph); })
    && std::all_of(phs.begin(), phs.end(), [&nary_boxed_chull] (const Poly& ph)
                   { return nary_boxed_chull.contains(ph); })
    && std::all_of(phs.begin(), phs.end(), [&seq_chull] (const Poly& ph)
                   { return seq_chull.contains(ph); });
  return ok;
}

namespace {

Poly
create_random_poly(dim_type dim, dim_type nr_cons) {
  Poly ph(dim);
  for (auto c = 0; c < nr_cons/2; ++c) {
    dim_type d1 = (rand()) % dim;
    dim_type d2 = (rand()) % dim;
    dim_type d3 = (rand()) % dim;
    dim_type c1 = rand() % 20;
    dim_type c2 = rand() % 20;
    dim_type c3 = rand() % 20;
    dim_type inhomo = rand() % 100;
    // Ensure to be polytope
    ph.add_con(c1 * Var(d1) + c2 * Var(d2) + c3 * Var(d3) >= -inhomo);
    ph.add_con(c1 * Var(d1) + c2 * Var(d2) + c3 * Var(d3) <= inhomo);
  }
  ph.minimize();
  return ph;
}

} // namespace

bool
test02() {
  dim_type nr = 10;
  dim_type nr_cons = 10;
  dim_type dim = 5;
  if (check_exp_eval()) {
    nr = 30;
    nr_cons = 30;
    dim = 15;
  }

  nout << "\nVector of " << nr << " random " << dim << "-dim polytopes, "
       << "each with " << nr_cons<< " cons." << std::endl;

  std::vector<Poly> phs;
  for (auto i = 0; i < nr; ++i) {
    phs.push_back(create_random_poly(dim, nr_cons));
  }

  Poly nary_chull(dim, Spec_Elem::EMPTY);
  {
    Clock clock;
    nary_chull = con_hull(phs.begin(), phs.end());
    nout << "  Time spent in n-ary con hull: ";
    clock.print_elapsed(nout);
    nout << std::endl;
  }

  Poly nary_boxed_chull(dim, Spec_Elem::EMPTY);
  {
    Clock clock;
    nary_boxed_chull = con_hull(phs.begin(), phs.end(), true);
    nout << "  Time spent in n-ary boxed con hull: ";
    clock.print_elapsed(nout);
    nout << std::endl;
  }

  Poly seq_chull(dim, Spec_Elem::EMPTY);
  {
    Clock clock;
    for (auto i = 0; i < nr; ++i) {
      seq_chull.con_hull_assign(phs[i]);
    }
    nout << "  Time spent in sequential con hull: ";
    clock.print_elapsed(nout);
    nout << std::endl;
  }

  bool ok = (seq_chull.contains(nary_chull)
             && nary_chull.contains(nary_boxed_chull));
  if (!ok)
    return false;

  if (!nary_chull.contains(seq_chull))
    nout << ">>> n-ary operator is more precise." << std::endl;

  if (!nary_boxed_chull.contains(nary_chull))
    nout << ">>> boxed operator is more precise." << std::endl;

  ok = std::all_of(phs.begin(), phs.end(), [&nary_chull] (const Poly& ph)
                   { return nary_chull.contains(ph); })
    && std::all_of(phs.begin(), phs.end(), [&nary_boxed_chull] (const Poly& ph)
                   { return nary_boxed_chull.contains(ph); })
    && std::all_of(phs.begin(), phs.end(), [&seq_chull] (const Poly& ph)
                   { return seq_chull.contains(ph); });

  return ok;
}

bool
test03() {
  dim_type nr = 10;
  dim_type nr_cons = 10;
  dim_type dim = 5;
  if (check_exp_eval()) {
    nr = 30;
    nr_cons = 30;
    dim = 20;
  }

  nout << "\nVector of " << nr << " replicas of the same "
       << dim << "-dim polytope "
       << "with " << nr_cons << " cons." << std::endl;

  Poly ph = create_random_poly(dim, nr_cons);
  std::vector<Poly> phs;
  for (auto i = 0; i < nr; ++i) {
    phs.push_back(ph);
  }

  Poly nary_chull(dim, Spec_Elem::EMPTY);
  {
    Clock clock;
    nary_chull = con_hull(phs.begin(), phs.end());
    nout << "  Time spent in n-ary con hull: ";
    clock.print_elapsed(nout);
    nout << std::endl;
  }

  Poly nary_boxed_chull(dim, Spec_Elem::EMPTY);
  {
    Clock clock;
    nary_boxed_chull = con_hull(phs.begin(), phs.end(), true);
    nout << "  Time spent in n-ary boxed con hull: ";
    clock.print_elapsed(nout);
    nout << std::endl;
  }

  Poly seq_chull(dim, Spec_Elem::EMPTY);
  {
    Clock clock;
    for (auto i = 0; i < nr; ++i) {
      seq_chull.con_hull_assign(phs[i]);
    }
    nout << "  Time spent in sequential con hull: ";
    clock.print_elapsed(nout);
    nout << std::endl;
  }

  // the n-ary con-hull is precise as the poly-hull,
  // if it is present in the set.
  bool ok = (nary_chull == ph
             && seq_chull == nary_chull
             && nary_chull == nary_boxed_chull);
  if (!ok)
    return false;

  ok = std::all_of(phs.begin(), phs.end(), [&nary_chull] (const Poly& ph)
                   { return nary_chull.contains(ph); })
    && std::all_of(phs.begin(), phs.end(), [&nary_boxed_chull] (const Poly& ph)
                   { return nary_boxed_chull.contains(ph); })
    && std::all_of(phs.begin(), phs.end(), [&seq_chull] (const Poly& ph)
                   { return seq_chull.contains(ph); });
  return ok;
}

bool
test04() {
  dim_type nr = 10;
  dim_type nr_cons = 10;
  dim_type dim = 5;
  if (check_exp_eval()) {
    nr = 30;
    nr_cons = 22;
    dim = 10;
  }

  std::vector<Poly> phs;

  nout << "\nVector of " << nr - 1 << " random "
            << dim << "-dim polytopes, "
            << "each with " << nr_cons << " cons, "
            << "plus their poly hull." << std::endl;

  for (auto i = 0; i < nr - 1; ++i) {
    phs.push_back(create_random_poly(dim, nr_cons));
  }
  // Add the poly_hull
  Poly hull(dim, Spec_Elem::EMPTY);
  {
    Clock clock;
    for (auto i = 0; i < nr - 1; ++i) {
      hull.poly_hull_assign(phs[i]);
    }
    hull.minimize();
    nout << "  Time spent in poly hull and minimize: ";
    clock.print_elapsed(nout);
    nout << std::endl;
  }
  phs.push_back(hull);

  Poly nary_chull(dim, Spec_Elem::EMPTY);
  {
    Clock clock;
    nary_chull = con_hull(phs.begin(), phs.end());
    nout << "  Time spent in n-ary con hull: ";
    clock.print_elapsed(nout);
    nout << std::endl;
  }

  Poly nary_boxed_chull(dim, Spec_Elem::EMPTY);
  {
    Clock clock;
    nary_boxed_chull = con_hull(phs.begin(), phs.end(), true);
    nout << "  Time spent in n-ary boxed con hull: ";
    clock.print_elapsed(nout);
    nout << std::endl;
  }

  Poly seq_chull(dim, Spec_Elem::EMPTY);
  {
    Clock clock;
    for (auto i = 0; i < nr; ++i) {
      seq_chull.con_hull_assign(phs[i]);
    }
    nout << "  Time spent in sequential con hull: ";
    clock.print_elapsed(nout);
    nout << std::endl;
  }

  // the n-ary con-hull is precise as the poly-hull,
  // if it is present in the set.
  bool ok = (nary_chull == hull
             && nary_chull == nary_boxed_chull
             && seq_chull.contains(nary_chull));

  if (!nary_chull.contains(seq_chull))
    nout << ">>> n-ary operator is more precise." << std::endl;

  return ok;
}

bool
test05() {

  Var A(0); Var B(1); Var C(2); Var D(3);
  Var E(4); Var F(5); Var G(6);

  Cons cs_x = { F == 0, D + E == 8, C == 8, A == 1, D - G >= 0,
                -4*B + 4*G >= 1, B >= 0, -G >= -1, -D >= -4 };
  Poly x(7);
  x.add_cons(cs_x);

  Cons cs_y = { D + E == 8, C == 8, B == 0, A == 1, -D >= -4, D >= 0,
                F >= 0, -F + G >= -4, -G >= -1, 4*G >= 1, -D - F + G >= -5 };
  Poly y(7);
  y.add_cons(cs_y);

  Poly kr = x;
  kr.con_hull_assign(y);

  print_cons(kr, "*** kr ***");

  con_hull(x, &y, &y + 1, false);

  print_cons(x, "*** x ***");

  bool ok = (x == kr);

  return ok;
}

// This test shows that con_hull happens to be more precise for U_Poly,
// due to an optimization that triggers its non-monotonicity.
bool test06() {
  dim_type dim = 3;

  Var A(0); Var B(1); Var C(2);

  // x_cs implies C >= 1
  Cons x_cs = {
    +A + 2*C >= 0,
    -A - C >= 1
  };

  // y_cs implies C >= 1
  Cons y_cs = {
    +B + 2*C >= 0,
    -B - C >= 1
  };

  Poly p_x(dim);
  p_x.add_cons(x_cs);
  p_x.minimize();
  U_Poly u_x(dim);
  u_x.add_cons(x_cs);
  u_x.minimize();

  Poly p_y(dim);
  p_y.add_cons(y_cs);
  p_y.minimize();
  U_Poly u_y(dim);
  u_y.add_cons(y_cs);
  u_y.minimize();

  auto p_res = p_x;
  p_res.con_hull_assign(p_y);
  p_res.minimize();

  auto u_res = u_x;
  u_res.con_hull_assign(u_y);
  u_res.minimize();

  print_cons(p_res, "*** After chull Poly X ***");
  print_cons(u_res, "*** After chull U_Poly X ***");

  Poly up_res(dim);
  up_res.add_cons(u_res.copy_cons());

  nout << std::endl;

  bool ok = true;
  if (p_res.contains(p_x) && p_res.contains(p_y))
    nout << "Poly result is an upper bound" << std::endl;
  else
    ok = false;
  if (u_res.contains(u_x) && u_res.contains(u_y))
    nout << "U_Poly result is an upper bound" << std::endl;
  else
    ok = false;

  if (p_res.contains(up_res))
    nout << "Poly result contains U_Poly result" << std::endl;
  else
    ok = false;
  if (not up_res.contains(p_res))
    nout << "U_Poly result does NOT contain Poly result" << std::endl;
  else
    ok = false;
  return ok;
}

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
  DO_TEST(test03);
  DO_TEST(test04);
  DO_TEST(test05);
  DO_TEST(test06);
END_MAIN
