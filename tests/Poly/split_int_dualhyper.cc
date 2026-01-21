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
#include "pplite_test.hh"
#include <vector>
#include <map>

namespace {

Gens
points_dual_hypercube(const dim_type dims,
                      const Linear_Expr& weight_center,
                      const Integer& half_diagonal) {
  Gens gs;
  for (dim_type axis = dims; axis-- > 0; ) {
    gs.push_back(point(weight_center + half_diagonal * Var(axis)));
    gs.push_back(point(weight_center - half_diagonal * Var(axis)));
  }
  return gs;
}

Poly
dual_hypercube(const dim_type dims,
               const Linear_Expr& weight_center,
               const Integer& half_diagonal) {
  Gens gs = points_dual_hypercube(dims, weight_center, half_diagonal);
  Poly ph(dims, Spec_Elem::EMPTY, Topol::CLOSED);
  ph.add_gens(gs);
  return ph;
}

void
build_polyhedra(const dim_type dims, std::vector<Poly>& ph) {
  Linear_Expr zero;
  Linear_Expr weight_center;

  // 1st-polyhedron.
  weight_center = zero;
  for (dim_type axis = dims; axis-- > 0; )
    weight_center += Var(axis);
  ph.push_back(dual_hypercube(dims, weight_center, 5));

  // 2nd-polyhedron.
  weight_center = zero;
  for (dim_type axis = dims; axis-- > 0; )
    weight_center += 2*Var(axis);
  ph.push_back(dual_hypercube(dims, weight_center, 4));

  // 3rd-polyhedron.
  weight_center = zero;
  for (dim_type axis = dims; axis-- > 0; )
    if (axis % 2 == 0)
      weight_center += 10*Var(axis);
    else
      weight_center += 2*Var(axis);
  ph.push_back(dual_hypercube(dims, weight_center, 5));

  // 4th-polyhedron.
  weight_center = zero;
  for (dim_type axis = dims; axis-- > 0; )
    if (axis % 2 == 0)
      weight_center += 10*Var(axis);
    else
      weight_center += Var(axis);
  ph.push_back(dual_hypercube(dims, weight_center, 4));
}

bool
test01() {
  std::vector<Poly> ph;

  dim_type first_dim = 4;
  dim_type last_dim = 5;

  bool ok = true;
  for (dim_type dims = first_dim; dims <= last_dim; dims++) {
    nout << endl
         << "++++++++ DIMENSIONS = " << dims << "  ++++++++"
         << endl;

    // Standard evaluation strategy.
    ph.clear();
    build_polyhedra(dims, ph);

    // For each poly in ph, check split with every constraint from another one.
    for (auto p : ph) {
      for (const auto& p1 : ph) {
        Cons cs = p1.copy_cons();
        for (auto c : cs) {
          Poly p_cache(p);

          Poly knres_then(p);
          Poly knres_else(p);
          auto [c_then, c_else] = detail::integral_complement_cons(c);
          knres_then.add_con(c_then);
          knres_else.add_con(c_else);

          Poly ph_else = p.integral_split(c);

          ok = (knres_then == p);
          if (!ok) {
            nout << "****** ph_cache ********* ";
            p_cache.ascii_dump(nout);
            nout << "****** ph_then ********* ";
            p.ascii_dump(nout);
            nout << "****** knres_then ********* ";
            knres_then.ascii_dump(nout);
            return ok;
          }
          ok = (knres_else == ph_else);
          if (!ok) {
            nout << "****** ph_cache ********* ";
            p_cache.ascii_dump(nout);
            nout << "****** ph_else ********* ";
            ph_else.ascii_dump(nout);
            nout << "****** knres_else ********* ";
            knres_else.ascii_dump(nout);
            return ok;
          }
        }
      }
    }
  }
  return ok;
}

bool
test02() {
  std::vector<Poly> ph;

  dim_type first_dim = 7;
  dim_type last_dim = 7;

  bool ok = true;
  for (dim_type dims = first_dim; dims <= last_dim; dims++) {
    nout << endl
         << "++++++++ DIMENSIONS = " << dims << "  ++++++++"
         << endl;

    // Standard evaluation strategy.
    ph.clear();
    build_polyhedra(dims, ph);

    Cons cs;
    cs.reserve(10);
    for (dim_type i = -4; i <= 5; ++i)
      cs.push_back(Var(0) <= i);

    for (const auto& poly : ph) {
      std::vector<Poly> ph_split;
      ph_split.resize(11, Poly(dims));
      std::vector<Poly> knres;
      knres.resize(11, Poly(dims));
      ph_split[0] = poly;
      knres[0] = poly;

      nout << "\nSplit mode \t | ";
      {
        Clock clock;
        for (dim_type i = 0; i < 10; ++i) {
          ph_split[i + 1] = ph_split[i].integral_split(cs[i]);
        }
        clock.print_elapsed(nout);
      }
      nout << "\nEager mode \t | ";
      {
        Clock clock;
        for (dim_type i = 0; i < 10; ++i) {
          auto [c_then, c_else] = detail::integral_complement_cons(cs[i]);
          knres[i + 1] = knres[i];
          knres[i].add_con(c_then);
          knres[i].minimize();
          knres[i + 1].add_con(c_else);
          knres[i + 1].minimize();
        }
        clock.print_elapsed(nout);
      }

      // Check
      for (dim_type i = 0; i <= 10; ++i) {
        ok = (ph_split[i] == knres[i]);
        if (!ok) {
          nout << " With constraint : ";
          cs[i].ascii_dump(nout);
          nout << "****** ph_then ********* ";
          ph_split[i].ascii_dump(nout);
          nout << "****** knres_then ********* ";
          knres[i].ascii_dump(nout);
          return ok;
        }
      }
      nout << "\n";
    }
  }
  return ok;
}

bool
test03() {
  // Build dual hypercube.
  auto build_dh = [](dim_type dims) {
    Linear_Expr zero;
    Linear_Expr weight_center = zero;
    for (dim_type axis = dims; axis-- > 0; )
      if (axis % 2 == 0)
        weight_center += 10*Var(axis);
      else
        weight_center += Var(axis);
    return dual_hypercube(dims, weight_center, 4);
  };

  dim_type first_dim = 8;
  dim_type last_dim = 8;
  if (check_exp_eval()) {
    first_dim = 8;
    last_dim = 11;
  }

  bool ok = true;
  for (dim_type dims = first_dim; dims <= last_dim; dims++) {

    nout << endl
         << "++++++++ DIMENSIONS = " << dims << "  ++++++++"
         << endl;

    Cons cs;
    cs.reserve(10);
    for (auto i = -4; i <= 5; ++i)
      cs.push_back(((i%2) ? Var(0) : Var(1)) < i);

    std::vector<Poly> phs(11, Poly(dims));
    phs[0] = build_dh(dims);

    std::vector<Poly> knres = phs;

    nout << "\nNew split: \t | ";
    {
      Clock clock;
      for (auto i = 0; i < 10; ++i)
        phs[i + 1] = phs[i].integral_split(cs[i]);
      clock.print_elapsed(nout);
    }
    nout << "\nOld split: \t | ";
    {
      Clock clock;
      for (auto i = 0; i < 10; ++i) {
        auto [c_then, c_else] = detail::integral_complement_cons(cs[i]);
        knres[i + 1] = knres[i];
        knres[i].add_con(c_then);
        knres[i].minimize();
        knres[i + 1].add_con(c_else);
        knres[i + 1].minimize();
      }
      clock.print_elapsed(nout);
    }

    // Check
    for (dim_type i = 0; i <= 10; ++i) {
      ok = (phs[i] == knres[i]);
      if (!ok) {
        nout << " With constraint : ";
        cs[i].ascii_dump(nout);
        nout << "****** ph_then ********* ";
        phs[i].ascii_dump(nout);
        nout << "****** knres_then ********* ";
        knres[i].ascii_dump(nout);
        return ok;
      }
    }
    nout << "\n";
  }

  return ok;
}

} // namespace

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
  DO_TEST(test03);
END_MAIN
