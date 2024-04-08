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
    gs.push_back(closure_point(weight_center - half_diagonal * Var(axis)));
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

long
computation(std::vector<Poly>& ph) {
  using IO_Operators::operator<<;
  const bool noisy
    = (getenv("PPLITE_NOISY_TESTS") != 0)
    || (getenv("PPLITE_VERY_NOISY_TESTS") != 0);

  nout << endl;
  nout << "working with 4 NNC dual hypercubes of dimension "
       << ph[0].space_dim() << endl;

  Clock clock;

  /**** Compute the intersection of ph[0] and ph[1]. ****/

  nout << "  - Computing intersection of ph[0] and ph[1]:" << endl;

  if (noisy) {
    // Print cardinalities of generator systems.
    Gens gs_0 = ph[0].copy_gens();
    Gens gs_1 = ph[1].copy_gens();
    nout << "    # ph[0].copy_gens() = " << gs_0.size() << endl;
    nout << "    # ph[1].copy_gens() = " << gs_1.size() << endl;
    vnout << "*** ph[0] generators ***" << endl
          << gs_0 << endl;
    vnout << "*** ph[1] generators ***" << endl
          << gs_1 << endl;

    // Print cardinalities of constraint systems.
    Cons cs_0 = ph[0].copy_cons();
    Cons cs_1 = ph[1].copy_cons();
    nout << "    # ph[0].copy_cons() = " << cs_0.size() << endl;
    nout << "    # ph[1].copy_cons() = " << cs_1.size() << endl;
    // Very noisy dump of arguments.
    vnout << "*** ph[0] constraints ***" << endl
          << cs_0 << endl;
    vnout << "*** ph[1] constraints ***" << endl
          << cs_1 << endl;
  }

  ph[0].intersection_assign(ph[1]);

  /**** Compute the intersection of ph[2] and ph[3]. ****/

  nout << "  - Computing intersection of ph[2] and ph[3]:" << endl;

  if (noisy) {
    // Print cardinalities of generator systems.
    Gens gs_2 = ph[2].copy_gens();
    Gens gs_3 = ph[3].copy_gens();
    nout << "    # ph[2].copy_gens() = " << gs_2.size() << endl;
    nout << "    # ph[3].copy_gens() = " << gs_3.size() << endl;
    // Very noisy dump of arguments.
    vnout << "*** ph[2] generators ***" << endl
          << gs_2 << endl;
    vnout << "*** ph[3] generators ***" << endl
          << gs_3 << endl;

    // Print cardinalities of constraint systems.
    Cons cs_2 = ph[2].copy_cons();
    Cons cs_3 = ph[3].copy_cons();
    nout << "    # ph[2].copy_cons() = " << cs_2.size() << endl;
    nout << "    # ph[3].copy_cons() = " << cs_3.size() << endl;
    // Very noisy dump of arguments.
    vnout << "*** ph[2] constraints ***" << endl
          << cs_2 << endl;
    vnout << "*** ph[3] constraints ***" << endl
          << cs_3 << endl;
  }

  ph[2].intersection_assign(ph[3]);

  /**** Compute the poly-hull of ph[0] and ph[2]. ****/
  nout << "  - Computing poly-hull of ph[0] and ph[2]:" << endl;

  if (noisy) {
    // Print cardinalities of generator systems.
    Gens gs_01 = ph[0].copy_gens();
    Gens gs_23 = ph[2].copy_gens();
    nout << "    # ph[0].copy_gens() = " << gs_01.size() << endl;
    nout << "    # ph[2].copy_gens() = " << gs_23.size() << endl;
    // Very noisy dump of arguments.
    vnout << "*** ph[0] generators ***" << endl
         << gs_01 << endl;
    vnout << "*** ph[2] generators ***" << endl
         << gs_23 << endl;
  }

  ph[0].poly_hull_assign(ph[2]);

  /**** Final conversion ****/

  Cons cs = ph[0].copy_cons();

  nout << "Final result timing: ";
  clock.print_elapsed(nout);
  nout << endl;

  // How many constraints obtained?
  const long cs_size = cs.size();

  // Print cardinality of final result.
  nout << "  - Final result is ph[0]:" << endl;
  nout << "    # ph[0].copy_cons() = " << cs_size << endl;
  // Very noisy dump of weakly-minimized final result.
  vnout << "*** ph[0] constraints ***" << endl;
  vnout << cs << endl;

  return cs_size;
}

bool
test01() {
  std::vector<Poly> ph;

  dim_type first_dim = 4;
  dim_type last_dim = 5;

  // Storing cardinalities of known results.
  using My_Map = std::map<std::pair<dim_type, int>, long>;
  My_Map::const_iterator known_result;
  My_Map cardinalities;

  using std::make_pair;

  cardinalities[make_pair(4, 25)] = 31;
  cardinalities[make_pair(4, 50)] = 41;
  cardinalities[make_pair(5, 25)] = 125;
  cardinalities[make_pair(5, 50)] = 150;

  int num_errors = 0;
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
      const long eval_card = computation(ph);

      // Check if there is a known result.
      known_result = cardinalities.find(make_pair(dims, percentage));
      if (known_result != cardinalities.end()
          && known_result->second != eval_card) {
        ++num_errors;
        nout << "Cardinality mismatch: "
             << "expected " << known_result->second << ", "
             << "obtained " << eval_card << ".\n";
      }
    }
  }
  return num_errors == 0;
}

} // namespace

BEGIN_MAIN
  DO_TEST(test01);
END_MAIN
