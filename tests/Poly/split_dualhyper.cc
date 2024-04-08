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
#include "pplite_test.hh"
#include <vector>
#include <map>

namespace {

void
closure_points_dual_hypercube(const dim_type dims,
                              const Linear_Expr& weight_center,
                              const Integer& half_diagonal,
                              Gens& gs) {
  // An ill-formed (it has no points at all) generator system
  // for a dual hypercube.
  for (dim_type axis = dims; axis-- > 0; ) {
    gs.push_back(closure_point(weight_center + half_diagonal * Var(axis)));
    // Let there be some ns
    gs.push_back(point(weight_center - half_diagonal * Var(axis)));
  }
}

void
add_facets(dim_type& to_be_added,
           Gens& gs,
           const Linear_Expr& expr,
           const dim_type axis,
           const dim_type dims,
           const Linear_Expr& weight_center,
           const Integer& half_diagonal) {
  // Return if we have already added all facets.
  if (to_be_added == 0)
    return;

  Linear_Expr expr1 = expr;
  expr1 += half_diagonal * Var(axis);
  Linear_Expr expr2 = expr;
  expr2 -= half_diagonal * Var(axis);

  if (axis == 0) {
    gs.push_back(point(dims * weight_center + expr1, dims));
    --to_be_added;
    if (to_be_added == 0)
      return;
    gs.push_back(point(dims * weight_center + expr2, dims));
    --to_be_added;
    return;
  }

  // Here axis > 0.
  // First recursive call with variable with index `axis'
  // having coordinate 1/dims.
  add_facets(to_be_added, gs, expr1,
             axis-1, dims, weight_center, half_diagonal);
  if (to_be_added == 0)
    return;
  // Second recursive call with variable with index `axis'
  // having coordinate -1/dims.
  add_facets(to_be_added, gs, expr2,
             axis-1, dims, weight_center, half_diagonal);
}

Poly
NNC_dual_hypercube(const dim_type dims,
                   const Linear_Expr& weight_center,
                   const Integer& half_diagonal,
                   const int facet_percentage) {
  Gens gs;
  closure_points_dual_hypercube(dims, weight_center, half_diagonal, gs);
  // Number of facets in the closed dual hypercube.
  dim_type num_facets = 1;
  for (dim_type axis = dims; axis-- > 0; )
    num_facets *= 2;
  dim_type facets_to_be_added = (num_facets * facet_percentage) / 100;
  if (facets_to_be_added == 0)
    // There has to be a point, at least.
    gs.push_back(point(weight_center));
  else
    add_facets(facets_to_be_added, gs, Linear_Expr(),
               dims-1, dims, weight_center, half_diagonal);
  // Actually build the polyhedron.
  Poly ph(dims, Spec_Elem::EMPTY, Topol::NNC);
  ph.add_gens(gs);
  return ph;
}

void
build_polyhedra(const dim_type dims,
                const int percentage,
                std::vector<Poly>& ph) {
  Linear_Expr zero;
  Linear_Expr weight_center;

  // 1st-polyhedron.
  weight_center = zero;
  for (dim_type axis = dims; axis-- > 0; )
    weight_center += Var(axis);
  ph.push_back(NNC_dual_hypercube(dims, weight_center, 5, percentage));

  // 2nd-polyhedron.
  weight_center = zero;
  for (dim_type axis = dims; axis-- > 0; )
    weight_center += 2*Var(axis);
  ph.push_back(NNC_dual_hypercube(dims, weight_center, 4, percentage));

  // 3rd-polyhedron.
  weight_center = zero;
  for (dim_type axis = dims; axis-- > 0; )
    if (axis % 2 == 0)
      weight_center += 10*Var(axis);
    else
      weight_center += 2*Var(axis);
  ph.push_back(NNC_dual_hypercube(dims, weight_center, 5, percentage));

  // 4th-polyhedron.
  weight_center = zero;
  for (dim_type axis = dims; axis-- > 0; )
    if (axis % 2 == 0)
      weight_center += 10*Var(axis);
    else
      weight_center += Var(axis);
  ph.push_back(NNC_dual_hypercube(dims, weight_center, 4, percentage));
}

Con
complement_con(const Con& beta) {
  auto expr = beta.linear_expr();
  auto inhomo = beta.inhomo_term();
  neg_assign(inhomo);
  neg_assign(expr);
  Con::Type type = Con::STRICT_INEQUALITY;
  if (beta.is_strict_inequality())
    type = Con::NONSTRICT_INEQUALITY;
  return Con(expr, inhomo, type);
}

bool
test01() {
  std::vector<Poly> ph;

  dim_type first_dim = 4;
  dim_type last_dim = 5;

  bool ok = true;
  for (dim_type dims = first_dim; dims <= last_dim; dims++) {
    for (int percentage = 25; percentage <= 50; percentage += 25) {

      nout << endl
           << "++++++++ DIMENSIONS = " << dims << "  ++++++++"
           << endl
           << "++++++++ PERCENTAGE = " << percentage << " ++++++++"
           << endl;

      // Standard evaluation strategy.
      ph.clear();
      build_polyhedra(dims, percentage, ph);

      // For each poly in ph, check split with every constraint from another one.
      for (auto p : ph) {
        for (const auto& p1 : ph) {
          Cons cs = p1.copy_cons();
          for (auto c : cs) {
            Poly p_cache(p);

            Poly knres_then(p);
            Poly knres_else(p);
            knres_then.add_con(c);
            knres_else.add_con(complement_con(c));

            Poly ph_else = p.split(c);

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
    for (int percentage = 25; percentage <= 50; percentage += 25) {

      nout << endl
           << "++++++++ DIMENSIONS = " << dims << "  ++++++++"
           << endl
           << "++++++++ PERCENTAGE = " << percentage << " ++++++++"
           << endl;

      // Standard evaluation strategy.
      ph.clear();
      build_polyhedra(dims, percentage, ph);

      Cons cs;
      cs.reserve(10);
      for (dim_type i = -4; i <= 5; ++i)
        cs.push_back(Var(0) <= i);

      for (const auto& poly : ph) {
        std::vector<Poly> ph_split;
        ph_split.resize(11, Poly(dims, Topol::NNC));
        std::vector<Poly> knres;
        knres.resize(11, Poly(dims, Topol::NNC));
        ph_split[0] = poly;
        knres[0] = poly;

        nout << "\nSplit mode \t | ";
        {
          Clock clock;
          for (dim_type i = 0; i < 10; ++i) {
            ph_split[i + 1] = ph_split[i].split(cs[i]);
          }
          clock.print_elapsed(nout);
        }
        nout << "\nEager mode \t | ";
        {
          Clock clock;
          for (dim_type i = 0; i < 10; ++i) {
            knres[i + 1] = knres[i];
            knres[i].add_con(cs[i]);
            knres[i].minimize();
            knres[i + 1].add_con(complement_con(cs[i]));
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
  }
  return ok;
}

bool
test03() {
  // Build dual hypercube.
  auto build_dh = [](dim_type dims, int percentage) {
    Linear_Expr zero;
    Linear_Expr weight_center = zero;
    for (dim_type axis = dims; axis-- > 0; )
      if (axis % 2 == 0)
        weight_center += 10*Var(axis);
      else
        weight_center += Var(axis);
    return (NNC_dual_hypercube(dims, weight_center, 4, percentage));
  };

  dim_type first_dim = 8;
  dim_type last_dim = 8;
  if (check_exp_eval()) {
    first_dim = 8;
    last_dim = 11;
  }

  bool ok = true;
  for (dim_type dims = first_dim; dims <= last_dim; dims++) {
    for (int percentage = 25; percentage <= 50; percentage += 25) {

      nout << endl
           << "++++++++ DIMENSIONS = " << dims << "  ++++++++"
           << endl
           << "++++++++ PERCENTAGE = " << percentage << " ++++++++"
           << endl;

      Cons cs;
      cs.reserve(10);
      for (auto i = -4; i <= 5; ++i)
          cs.push_back(((i%2) ? Var(0) : Var(1)) < i);

      std::vector<Poly> phs(11, Poly(dims, Topol::NNC));
      phs[0] = build_dh(dims, percentage);

      std::vector<Poly> knres = phs;

      nout << "\nNew split: \t | ";
      {
        Clock clock;
        for (auto i = 0; i < 10; ++i)
          phs[i + 1] = phs[i].split(cs[i]);
        clock.print_elapsed(nout);
      }
      nout << "\nOld split: \t | ";
      {
        Clock clock;
        for (auto i = 0; i < 10; ++i) {
          knres[i + 1] = knres[i];
          knres[i].add_con(cs[i]);
          knres[i].minimize();
          knres[i + 1].add_con(complement_con(cs[i]));
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
  }

  return ok;
}

} // namespace

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
  DO_TEST(test03);
END_MAIN
