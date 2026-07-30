/* PPLite: a lightweight library for convex polyhedra derived from PPL.
   Copyright (C) 2026 Enea Zaffanella <enea.zaffanella@unipr.it>

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

#include "pplite-config.h"
#include "utils.hh"
#include "clock.hh"
#include "Bits.hh"
#include "Gen.hh"
#include "Integer.hh"
#include "Rational.hh"
#include "Poly.hh"
#include "Poly_min.hh"
#include "Poly_templ.hh"
#include "Var.hh"

#include <algorithm>
#include <iostream>
#include <memory>
#include <optional>
#include <set>
#include <stack>
#include <vector>

namespace pplite {

namespace minkowski {

/*
  Returns adjacency info for the points of (closed, non-empty) polytope ph.
  Let adj be the vector returned: for i and j indices of points
  in the skeleton generator system of ph, which is assumed
  to be minimized, adj[i][j] is set iff i and j are adjacent.
*/
std::vector<Bits>
get_adj_vertices(const Poly_Impl& ph) {
  assert(ph.space_dim() > 0 && not ph.is_empty());
  assert(ph.is_topologically_closed() && ph.is_bounded());
  assert(ph.is_minimized());
  const auto& gs = ph.gs.sk_rows;
  const auto& sat = ph.sat_c;
  const dim_type min_sat = ph.space_dim() - 1;
  const dim_type num_eqs = num_rows(ph.cs.sing_rows);
  const dim_type max_sat = sat.num_cols() + num_eqs;
  auto sz = num_rows(gs);
  std::vector<Bits> adj(sz);
  for (auto i = 0; i != sz; ++i) {
    assert(gs[i].is_point());
    const dim_type i_num_ones = sat[i].count_ones();
    for (auto j = i + 1; j != sz; ++j) {
      assert(gs[j].is_point());
      // quick non-adj test
      const dim_type new_satrow_ones = sat[i].count_ones_in_union(sat[j]);
      if (detail::quick_non_adj_test(min_sat, max_sat, new_satrow_ones))
        continue;
      // quick adj test
      const auto j_num_ones = sat[j].count_ones();
      const auto max_ones = std::max(i_num_ones, j_num_ones);
      if (max_ones + 1 == new_satrow_ones) {
        // adjacent
      } else {
        // (slow) combinatorial test
        Bits new_satrow = Bits::get_union(sat[i], sat[j]);
        if (not detail::combinatorial_adj_test(0, sz, i, j, new_satrow, sat))
          continue;
      }
      // i and j are adjacent points
      adj[i].set(j);
      adj[j].set(i);
    }
  }
  return adj;
}

Gens
get_points(const Poly_Impl& ph) {
  Gens res;
  for (const auto& g : ph.gens()) {
    if (g.is_point())
      res.push_back(g);
  }
  return res;
}

Gens
get_lines_and_rays(const Poly_Impl& ph) {
  Gens res;
  for (const auto& g : ph.gens()) {
    if (g.is_line_or_ray())
      res.push_back(g);
  }
  return res;
}

Poly_Impl
get_polytope(const Poly_Impl& x) {
  assert(not x.is_empty());
  assert(x.is_topologically_closed() && x.is_minimized());
  auto res = Poly_Impl(x.space_dim(), Spec_Elem::EMPTY, x.topology());
  res.add_gens(get_points(x));
  res.minimize();
  return res;
}

Gen
get_edge_direction(const Gen& px, const Gen& py) {
  assert(px.is_point() && py.is_point());
  assert(px != py);
  return ray(py.linear_expr() * px.divisor()
             - px.linear_expr() * py.divisor());
}

Gen sum_points(const Gen& px, const Gen& py) {
  assert(px.is_point() && py.is_point());
  return point(px.linear_expr() * py.divisor()
               + py.linear_expr() * px.divisor(),
               px.divisor() * py.divisor());
}

// a vertex in Minkowski's sum is "composed" by adding
// the i-th vertex in 1st arg with the j-th vertex in 2nd arg;
// hence, we represent it using indices (i,j)
using VComp = std::pair<dim_type, dim_type>;

Gen VComp2Gen(const VComp& v, const Gens& gs1, const Gens& gs2) {
  const auto& [i, j] = v;
  assert(0 <= i && i < num_rows(gs1));
  assert(0 <= j && j < num_rows(gs2));
  return sum_points(gs1[i], gs2[j]);
}

// The root is the vertex maximizing any (generic!) linear function.
VComp
get_root_composition(const Gens& gs1, const Gens& gs2) {

  // we arbitrarily select to maximize expression
  //    Var(0) + Var(1) + ... + Var(sdim-1)
  // i.e., just sum coefficients, no need to compute scalar product
  auto get_value = [](const Gen& g) {
    Integer numer;
    for (const auto& coeff : g.linear_expr())
      numer += coeff;
    return Rational(numer, g.divisor());
  };

  auto maximizing_point = [&get_value](const Gens& gs) {
    auto it = gs.begin();
    auto it_end = gs.end();
    assert(it != it_end && it->is_point());

    auto max_it = it;
    auto max_val = get_value(*it);
    for (++it; it != it_end; ++it) {
      assert(it->is_point());
      auto val = get_value(*it);
      if (val > max_val) {
        max_it = it;
        max_val = std::move(val);
      }
    }
    return std::distance(gs.begin(), max_it);
  };

  return VComp { maximizing_point(gs1), maximizing_point(gs2) };
}

// Returns the sum of (closed, non-empty and minimized) *polytopes* x and y
Poly_Impl
polytope_sum(const Poly_Impl& x, const Poly_Impl& y) {
  assert(x.space_dim() > 0);
  assert(x.space_dim() == y.space_dim());
  assert(x.is_topologically_closed() && y.is_topologically_closed());
  assert(not x.is_empty() && not y.is_empty());
  assert(x.is_minimized() && y.is_minimized());
  // x and y are polytopes
  assert(x.is_bounded() && y.is_bounded());

  // we only work on the skel component
  const auto& x_gs = x.gs.sk_rows;
  const auto& y_gs = y.gs.sk_rows;

  const auto x_adj = get_adj_vertices(x);
  const auto y_adj = get_adj_vertices(y);

  // comparison function for associative container
  auto gen_cmp = [](const Gen& g1, const Gen& g2) {
    return compare(g1, g2) < 0;
  };
  // use a set to efficiently check for duplicates
  std::set<Gen, decltype(gen_cmp)> sum_vertices(gen_cmp);

  auto is_duplicate_vertex = [&sum_vertices](const Gen& g) {
    return sum_vertices.find(g) != sum_vertices.end();
  };

  VComp root_v = get_root_composition(x_gs, y_gs);
  sum_vertices.insert(VComp2Gen(root_v, x_gs, y_gs));

  std::stack<VComp> working;
  working.push(root_v);

  /*
    Some adj_oracle checks are "symmetric"; if any call
        adj_oracle( (xi, yi), true, adj1)
    tests the adjacency of next_v = (adj1, adj2) to curr_v = (xi, yi)
    and we have *both* xi != adj1 and yi != adj2,
    then the symmetric call
        adj_oracle( (xi, yi), false, adj2),
    will test the adjacency of the same next_v
    (similarly if the order of the oracle calls is reversed);
    in both cases, the second oracle test is useless,
    either because next_v is not adjacent to curr_v or because
    it would anyway produce a duplicate vertex of the minkowski sum.
    To avoid wasting time redoing these symmetric tests,
    we keep track (in a set) of the symmetric tests already done.

    NOTE: this cannot be replaced by the check is_duplicate_vertex,
    because the latter can only succeed when next_v *is* adjacent to
    curr_v; the symmetric test is able to save repeated oracle checks
    returning nullopt (i.e., non-adjacent).
  */
  using Symm_Oracle_Test = std::tuple<VComp, bool, dim_type>;
  using Symm_Oracle_Tests = std::set<Symm_Oracle_Test>;
  Symm_Oracle_Tests symm_oracle_tests;

  auto record_symm_oracle_test
    = [&symm_oracle_tests](const VComp& curr_v, bool from_x, dim_type adj1) {
      symm_oracle_tests.insert(Symm_Oracle_Test{ curr_v, from_x, adj1 });
    };

  auto is_duplicate_oracle_test
    = [&symm_oracle_tests](const VComp& curr_v, bool from_x, dim_type adj2) {
      Symm_Oracle_Test t { curr_v, from_x, adj2 };
      return symm_oracle_tests.find(t) != symm_oracle_tests.end();
    };

  /*
    Note: this adjacency oracle for Minkowski sum is (loosely) based on
    a paper by K. Fukuda "From the zonotope construction to the Minkowski
    addition of convex polytopes", J. Symb. Comput., 2004
    (from_x and adj1 can be seen to encode the Delta-counter)
  */
  auto adj_oracle = [&](const VComp& curr_v,
                        bool from_x, dim_type adj1
                        ) -> std::optional<VComp> {
    const auto& [xi, yi] = curr_v;
    /* Note:
       from_x = true means that (xi, adj1) is a pair of indices
       of adjacent points indexing the x_gs component;
       from_x = false means that (yi, adj1) is a pair of indices
       of adjacent points indexing the y_gs component
    */
    assert((from_x && xi != adj1 && x_adj[xi][adj1])
           ||
           (not from_x && yi != adj1 && y_adj[yi][adj1]));

    // use subscript 1 or 2 for the x/y components
    auto idx1 = from_x ? xi : yi;
    auto idx2 = from_x ? yi : xi;
    const auto& gs1 = from_x ? x_gs : y_gs;
    const auto& gs2 = from_x ? y_gs : x_gs;
    const auto& adj_mat1 = from_x ? x_adj : y_adj;
    const auto& adj_mat2 = from_x ? y_adj : x_adj;

    // compute edge from v1
    auto edge1 = get_edge_direction(gs1[idx1], gs1[adj1]);

    Gens nonpar_edges;
    // collect all edges (starting from v2) *not* parallel to edge1;
    // also track the index of a parallel edge, if found
    dim_type v2_par_index = -1;
    for (auto adj2 : adj_mat2[idx2]) {
      auto edge2 = get_edge_direction(gs2[idx2], gs2[adj2]);
      // since we represent edges using (normalized) rays,
      // they are parallel only when they are equal
      if (edge1 != edge2)
        nonpar_edges.push_back(std::move(edge2));
      else {
        v2_par_index = adj2;
        // Here both idx1 != adj1 and idx2 != adj2, meaning that
        // this oracle test can be made symmetrically:
        // check if the symmetric test was already computed.
        if (is_duplicate_oracle_test(curr_v, not from_x, adj2))
          return std::nullopt;
        else
          record_symm_oracle_test(curr_v, from_x, adj1);
      }
    }

    // if the adjacency test is successful, we will return next_v
    auto adj2 = (v2_par_index >= 0) ? v2_par_index : idx2;
    VComp next_v = from_x ? VComp{adj1, adj2} : VComp{adj2, adj1};
    // even before completing the test for adjacency (which is expensive)
    // do check if next_v would generate a vertex already seen.
    if (is_duplicate_vertex(VComp2Gen(next_v, x_gs, y_gs)))
      return std::nullopt;

    // now add to nonpar_edges all edges starting from v1 (except edge1)
    for (auto other_adj1 : adj_mat1[idx1]) {
      if (other_adj1 == adj1)
        continue;
      auto other_edge1 = get_edge_direction(gs1[idx1], gs1[other_adj1]);
      nonpar_edges.push_back(std::move(other_edge1));
    }

    // remove possible duplicates (worth it?)
    std::sort(nonpar_edges.begin(), nonpar_edges.end(), gen_cmp);
    nonpar_edges.erase(std::unique(nonpar_edges.begin(), nonpar_edges.end()),
                       nonpar_edges.end());

    // Here we could (should!) check feasibility of an LP;
    // rather, we minimize by (inefficiently) checking all generators
    auto lp = Poly_Impl(x.space_dim(), Spec_Elem::UNIVERSE, Topol::CLOSED);
    for (const auto& npe : nonpar_edges)
      lp.add_con(npe.linear_expr() >= 0);
    const auto& expr1 = edge1.linear_expr();
    lp.add_con(expr1 >= -1); // arbitrary negative value
    lp.add_con(expr1 <= 0);

    Rational value;
    bool has_min = lp.min(Affine_Expr(expr1), value);
    bool is_adj = has_min && sgn(value) == -1;
    if (is_adj)
      return next_v;
    else
      return std::nullopt;
  }; // lambda adj_oracle

  auto check_adj = [&](const VComp& curr_v, bool from_x, dim_type adj_idx) {
    auto opt_next_v = adj_oracle(curr_v, from_x, adj_idx);
    if (opt_next_v) {
      auto next_v = opt_next_v.value();
      working.push(next_v);
      auto g = VComp2Gen(next_v, x_gs, y_gs);
      assert(not is_duplicate_vertex(g));
      sum_vertices.insert(std::move(g));
    }
  };

  while (not working.empty()) {
    auto curr_v = working.top();
    working.pop();
    auto [xi, yi] = curr_v;
    // case 1: check edges starting from xi (from_x = true)
    for (auto adj_xi : x_adj[xi])
      check_adj(curr_v, true, adj_xi);
    // case 2: check edges starting from yi (from_x = false)
    for (auto adj_yi : y_adj[yi])
      check_adj(curr_v, false, adj_yi);
  }

  auto res = Poly_Impl(x.space_dim(), Spec_Elem::EMPTY, x.topology());
  res.add_gens(sum_vertices.begin(), sum_vertices.end());

#ifndef NDEBUG
  // the algorithm should not produce redundant vertices
  dim_type card_adding = sum_vertices.size();
  auto res_copy = res; // working on a copy (it triggers minimization)
  dim_type card_res = res_copy.num_min_gens();
  assert(card_adding == card_res &&
         "ERROR in Minkowski's sum: adjacency oracle is returning "
         "a redundant set of vertices");
#endif

  return res;
}

/* only for debugging purposes: computationally very heavy */
Poly_Impl
naif_polytope_sum(const Poly_Impl& x, const Poly_Impl& y) {
  assert(x.space_dim() > 0);
  assert(x.space_dim() == y.space_dim());
  assert(x.is_topologically_closed() && y.is_topologically_closed());
  assert(not x.is_empty() && not y.is_empty());
  assert(x.is_minimized() && y.is_minimized());
  // x and y are polytopes
  assert(x.is_bounded() && y.is_bounded());

  Gens sum_gs;
  for (const auto& x_p : x.gens()) {
    assert(x_p.is_point());
    for (const auto& y_p : y.gens()) {
      assert(y_p.is_point());
      sum_gs.push_back(sum_points(x_p, y_p));
    }
  }
  auto res = Poly_Impl(x.space_dim(), Spec_Elem::EMPTY, x.topology());
  res.add_gens(sum_gs);
  res.minimize(); // not really needed
  return res;
}

} // namespace minkowski

Poly_Impl
Poly_Impl::minkowski_sum(const Poly_Impl& y) const {
  const auto& x = *this;
  assert(x.space_dim() == y.space_dim());
  // CHECK ME: is it worth having a version working on polyhedra
  // that are not topologically closed? Probably, no.
  assert(x.is_topologically_closed() && y.is_topologically_closed());

  auto res = Poly_Impl(x.space_dim(), Spec_Elem::EMPTY, x.topology());

  if (x.is_empty() || y.is_empty())
    return res;
  if (x.is_universe() || y.is_universe()) {
    res.set_universe();
    return res;
  }

  assert(x.space_dim() > 0);
  x.minimize();
  y.minimize();

  // Note: here we could check if any of x or y is a singleton point,
  // in which case the sum becomes an affine translation; we disregard
  // this potential optimization, as it is an unlikely case.

  if (x.is_bounded() && y.is_bounded()) {
    // x and y are polytopes
    res = minkowski::polytope_sum(x, y);
#if 0 // only for debugging
    auto res_naif = minkowski::naif_polytope_sum(x, y);
    assert(res.equals(res_naif));
#endif
    return res;
  }

  // x or y is unbounded (i.e., has lines/rays):
  // take a copy and add other argument lines/rays to it
  // (this is done eagerly so as to possibly reduce number of points)
  auto xx = x;
  xx.add_gens(minkowski::get_lines_and_rays(y));
  xx.minimize();
  if (xx.is_universe())
    return xx;
  auto yy = y;
  yy.add_gens(minkowski::get_lines_and_rays(x));
  yy.minimize();
  if (yy.is_universe())
    return yy;

  // extract polytopes from xx and yy
  auto x_tope = minkowski::get_polytope(xx);
  auto y_tope = minkowski::get_polytope(yy);
  res = minkowski::polytope_sum(x_tope, y_tope);
#if 0 // only for debugging
  auto res_naif = minkowski::naif_polytope_sum(x_tope, y_tope);
  assert(res.equals(res_naif));
#endif

  // add lines/rays (note: xx already includes those in yy)
  res.add_gens(minkowski::get_lines_and_rays(xx));
  return res;
}

} // namespace pplite
