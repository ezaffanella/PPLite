/* PPLite: a lightweight library for convex polyhedra derived from PPL.
   Copyright (C) 2001-2010 Roberto Bagnara <bagnara@cs.unipr.it>
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

#include "pplite-config.h"
#include "utils.hh"
#include "Poly.hh"
#include "Bits.hh"
#include "Integer.hh"
#include "Scalar_Prod.hh"
#include "Poly_templ.hh"
#include "Poly_min.hh"

namespace pplite {

bool
Poly_Impl::join_assign_if_exact(const Poly_Impl& y) {
  auto& x = *this;
  assert(x.topology() == y.topology());
  assert(x.space_dim() == y.space_dim());
  return (x.is_topologically_closed() && y.is_topologically_closed())
    ? x.closed_join_assign_if_exact(y)
    : x.nnc_join_assign_if_exact(y);
}

namespace {

Bits
subsumed_by(const Gens& gs, const Poly_Impl& ph) {
  Bits res;
  for (auto i : bwd_index_range(gs)) {
    if (ph.relation_with(gs[i]).implies(Poly_Gen_Rel::subsumes()))
      res.set(i);
  }
  return res;
}

Bits
non_redundant_in(const Gens& gs, const Poly_Impl& ph) {
  Bits res = subsumed_by(gs, ph);
  res.complement_until(num_rows(gs));
  return res;
}

Bits
points_non_redundant_in_closure(const Gens& gs, const Poly_Impl& ph) {
  Bits res;
  if (ph.is_topologically_closed()) {
    for (auto i : bwd_index_range(gs)) {
      const auto& g = gs[i];
      if (not g.is_point())
        continue;
      if (not ph.relation_with(g).implies(Poly_Gen_Rel::subsumes()))
        res.set(i);
    }
  } else {
    // ph not topologically closed, we cannot use method relation_with
    const auto& cs = ph.cons();
    for (auto i : bwd_index_range(gs)) {
      const auto& g = gs[i];
      if (not g.is_point())
        continue;
      for (const auto& c : cs) {
        auto sp_sign = sp::sign(c, g);
        if (sp_sign < 0 || (sp_sign > 0 && c.is_equality())) {
          res.set(i);
          break;
        }
      }
    }
  }
  return res;
}

} // namespace

/*
  Algorithm based on Theorem 3.2 in BHZ CGTA 2010:
  ph1 hull ph2 is *not* exact iff
  there exist c1 in cs1 and g1 in gs1 s.t.
    (A) g1 saturates c1,
    (B) ph2 violates c1, and
    (C) ph2 does not subsume g1.
  For efficiency, implementation also performs incomplete tests
  using bounding boxes and affine dimensions.
*/
bool
Poly_Impl::closed_join_assign_if_exact(const Poly_Impl& y) {
  Poly_Impl& x = *this;
  // Check preconditions
  assert(x.is_topologically_closed());
  assert(y.is_topologically_closed());
  assert(x.space_dim() == y.space_dim());

  // 0-dim case
  if (x.space_dim() == 0) {
    x.poly_hull_assign(y);
    return true;
  }
  // empty argument(s) case
  y.minimize();
  if (y.is_empty())
    return true;
  x.minimize();
  if (x.is_empty()) {
    x = y;
    return true;
  }

  // Disjoint bounding boxes case: hull not exact.
  // This is a quick (correct, but incomplete) check.
  BBox x_box = x.get_bounding_box();
  BBox y_box = y.get_bounding_box();
  if (x_box.is_disjoint_from(y_box))
    return false;

  // If affine dimensions differ, then join is exact
  // only when inclusion holds (note: correct for closed polyhedra).
  if (x.affine_dim() < y.affine_dim()) {
    if (not y_box.contains(x_box))
      return false;
    if (not y.boxed_contains(x))
      return false;
    x = y;
    return true;
  } else if (x.affine_dim() > y.affine_dim()) {
    if (not x_box.contains(y_box))
      return false;
    return x.boxed_contains(y);
  }

  assert(x.affine_dim() == y.affine_dim());

  // Complexity bound in Theorem 3.2 is asymmetric:
  // let ph1 be the "smaller" argument, where the metric
  // is the product of cons/gens cardinalities.
  auto x_sz = x.num_min_cons() * x.num_min_gens();
  auto y_sz = y.num_min_cons() * y.num_min_gens();
  bool x_is_ph1 = (x_sz <= y_sz);
  const auto& ph1 = x_is_ph1 ? x : y;
  const auto& ph2 = x_is_ph1 ? y : x;
  // Already mimimized, no pending cs/gs
  const auto& cs1 = ph1.cs;
  const auto& gs1 = ph1.gs;
  const auto& gs2 = ph2.gs;

  auto ln_sub1 = subsumed_by(gs1.sing_rows, ph2);
  auto sk_sub1 = subsumed_by(gs1.sk_rows, ph2);
  auto num_sub1 = ln_sub1.size() + sk_sub1.size();
  auto num_gs1 = num_rows(gs1.sing_rows) + num_rows(gs1.sk_rows);
  if (num_sub1 == num_gs1) {
    // ph1 included in ph2: ph2 is exact join
    // (no g1 fullfilling condition C)
    if (x_is_ph1)
      x = y;
    return true;
  }

  auto ln_sub2 = subsumed_by(gs2.sing_rows, ph1);
  auto sk_sub2 = subsumed_by(gs2.sk_rows, ph1);
  auto num_sub2 = ln_sub2.size() + sk_sub2.size();
  auto num_gs2 = num_rows(gs2.sing_rows) + num_rows(gs2.sk_rows);
  if (num_sub2 == num_gs2) {
    // ph2 included in ph1: ph1 is exact join
    // (no c1 fullfilling condition B)
    if (not x_is_ph1)
      x = y;
    return true;
  }

  // Here we know that no inclusion relation holds,
  // hence properties A and B are met by some c1 and g1;
  // if sub1 is empty, then property C holds for all of gs1,
  // so that the hull is not exact;
  // same for sub2 (exchanging roles of ph1 and ph2).
  if (num_sub1 == 0 || num_sub2 == 0)
    return false;

  // End of tricks: have to check properties A, B, C.
  // We start by checking equalities (in cs1.sing_rows):
  // in this case, saturation condition A always holds;
  if (num_rows(cs1.sing_rows) > 0) {
    // The checks for B and C are independent;
    // first check for property C
    bool C_holds
      = (ln_sub1.size() < num_rows(gs1.sing_rows))
      || (sk_sub1.size() < num_rows(gs1.sk_rows));
    // Then (if needed) check equalities for property B
    bool B_holds = false;
    if (C_holds) {
      for (const auto& c1 : cs1.sing_rows) {
        if (not ph2.relation_with(c1)
            .implies(Poly_Con_Rel::is_included())) {
          B_holds = true;
          break;
        }
      }
    }
    if (B_holds && C_holds)
      return false;
  }

  // Now we check inequalities (in cs1.sk_rows).
  if (num_rows(cs1.sk_rows) > 0) {
    // check C on singular generators (lines).
    bool C_holds_for_a_line = (ln_sub1.size() < num_rows(gs1.sing_rows));
    // For inequalities, saturation info is stored in ph1.sat_g;
    // the checks are not fully independent.
    const auto& sat1 = ph1.sat_g;
    assert(sat1.num_rows() == num_rows(cs1.sk_rows));
    const auto num_sk_gs1 = num_rows(gs1.sk_rows);
    for (auto i : bwd_index_range(cs1.sk_rows)) {
      const auto& c1 = cs1.sk_rows[i];
      if (not ph2.relation_with(c1)
          .implies(Poly_Con_Rel::is_included())) {
        // B holds.
        // If C was holding for a line, we are done
        // (lines are always saturated, hence A holds too).
        if (C_holds_for_a_line)
          return false;
        // Check for A and C on skel gens.
        // gens in sat1[i] do NOT saturate c1 (A does not hold);
        // gens in sk_sub1 are subsumed by ph2 (C does not hold);
        // if their union does NOT cover all of gs1.sk_rows,
        // then we found g1 satisfying both A and C.
        auto count = sk_sub1.count_ones_in_union(sat1[i]);
        if (count < num_sk_gs1)
          return false;
      }
    }
  }

  // The join is exact: compute it using the knowledge
  // that some of the gens are redundant
  // (and all gens in ln_sub1/ln_sub2 are redundant).
  const auto& y_gs = x_is_ph1 ? gs2.sk_rows : gs1.sk_rows;
  const auto& y_sub = x_is_ph1 ? sk_sub2 : sk_sub1;
  for (auto i : index_range(y_gs)) {
    if (not y_sub[i])
      x.add_gen(y_gs[i]);
  }
  return true;
}

/*
  Algorithm based on Theorem 3.6 in BHZ CGTA 2010:
  ph1 hull ph2 is *not* exact iff for i,j in { 1, 2 }, where i != j,
  there exist c_i in cs_i and g_i in gs_i s.t.
    (A) g_i saturates c_i,
    (B) ph_j violates c_i,
    and any one of the following holds:
    (C1) g_i is a ray or closure point not subsumed by ph_j;
      or
    (C2) g_i is a point, c_i is non-strict and g_i is not subsumed
         by the topological closure of ph_j;
      or
    (C3) c_i is strict and saturated by point p in (ph_i \phull ph_j)
         which is not subsumed by p_j.
*/
bool
Poly_Impl::nnc_join_assign_if_exact(const Poly_Impl& y) {
  Poly_Impl& x = *this;
  // Check preconditions
  assert(x.topology() == Topol::NNC);
  assert(y.topology() == Topol::NNC);
  assert(x.space_dim() == y.space_dim());

  // 0-dim case
  if (x.space_dim() == 0) {
    x.poly_hull_assign(y);
    return true;
  }
  // empty argument(s) case
  y.minimize();
  if (y.is_empty())
    return true;
  x.minimize();
  if (x.is_empty()) {
    x = y;
    return true;
  }

  // Disjoint bounding boxes case: hull not exact.
  // This is a quick (correct, but incomplete) check.
  // It is correct because boxes are topologically closed.
  BBox x_box = x.get_bounding_box();
  BBox y_box = y.get_bounding_box();
  if (x_box.is_disjoint_from(y_box))
    return false;

  // Everything is mimimized (no pending cs/gs)

  // materialize x nonskel points
  Gens x_gs_ns;
  x_gs_ns.reserve(num_rows(x.gs.ns_rows));
  for (const auto& ns : x.gs.ns_rows) {
    x_gs_ns.push_back(detail::materialize(ns, x.gs.sk_rows));
  }

  // compute non-redundant bitsets for x
  auto x_ln_nonred_in_y = non_redundant_in(x.gs.sing_rows, y);
  auto x_sk_nonred_in_y = non_redundant_in(x.gs.sk_rows, y);
  auto x_ns_nonred_in_y = non_redundant_in(x_gs_ns, y);

  // test for inclusion
  if (x_ln_nonred_in_y.empty() &&
      x_sk_nonred_in_y.empty() &&
      x_ns_nonred_in_y.empty()) {
    // x is included in y: exact hull is y
    x = y;
    return true;
  }

  // materialize y nonskel points
  Gens y_gs_ns;
  y_gs_ns.reserve(num_rows(y.gs.ns_rows));
  for (const auto& ns : y.gs.ns_rows) {
    y_gs_ns.push_back(detail::materialize(ns, y.gs.sk_rows));
  }

  // compute non-redundant bitsets for y
  auto y_ln_nonred_in_x = non_redundant_in(y.gs.sing_rows, x);
  auto y_sk_nonred_in_x = non_redundant_in(y.gs.sk_rows, x);
  auto y_ns_nonred_in_x = non_redundant_in(y_gs_ns, x);

  // test for inclusion
  if (y_ln_nonred_in_x.empty() &&
      y_sk_nonred_in_x.empty() &&
      y_ns_nonred_in_x.empty()) {
    // y is included in x: exact hull is x
    return true;
  }

  // further divide x non-red gens between points and non-points;
  // points should be checked wrt topological closure of y
  Bits x_sk_np_nonred_in_y;
  for (auto i : x_sk_nonred_in_y) {
    if (not x.gs.sk_rows[i].is_point())
      x_sk_np_nonred_in_y.set(i);
  }
  Bits x_sk_p_nonred_in_y_closure
    = points_non_redundant_in_closure(x.gs.sk_rows, y);
  Bits x_ns_nonred_in_y_closure
    = points_non_redundant_in_closure(x_gs_ns, y);

  // same for y
  Bits y_sk_np_nonred_in_x;
  for (auto i : y_sk_nonred_in_x) {
    if (not y.gs.sk_rows[i].is_point())
      y_sk_np_nonred_in_x.set(i);
  }
  Bits y_sk_p_nonred_in_x_closure
    = points_non_redundant_in_closure(y.gs.sk_rows, x);
  Bits y_ns_nonred_in_x_closure
    = points_non_redundant_in_closure(y_gs_ns, x);


  // checks conditions C1 and C2 on x equality constraints;
  // returns *true* if union is *NOT* exact, false if test not conclusive.
  auto check_equalities
    = [](const Poly_Impl& x, const Poly_Impl& y,
         const Bits& x_ln_nonred_in_y,
         const Bits& x_sk_np_nonred_in_y,
         const Bits& x_sk_p_nonred_in_y_closure,
         const Bits& x_ns_nonred_in_y_closure) {
    bool C1_holds
      = not x_ln_nonred_in_y.empty()
      || not x_sk_np_nonred_in_y.empty();
    bool C2_holds
      = not x_sk_p_nonred_in_y_closure.empty()
      || not x_ns_nonred_in_y_closure.empty();
    // No need to check for C3 (not a strict inequality)
    if (not C1_holds && not C2_holds)
      return false;
    // check for B, i.e., an equality in x violating y;
    // the equality saturates all gens of x, hence A will also hold.
    for (const auto& c : x.cs.sing_rows) {
      if (not y.relation_with(c).implies(Poly_Con_Rel::is_included()))
        return true;
    }
    return false;
  };

  // checks conditions C1 and C2 (but not C3) on x skel constraints;
  // returns true if union is NOT exact; false if test inconclusive.
  auto check_C1_C2_on_skel_ineqs
    = [](const Poly_Impl& x, const Poly_Impl& y,
         const Bits& x_ln_nonred_in_y,
         const Bits& x_sk_np_nonred_in_y,
         const Bits& x_sk_p_nonred_in_y_closure,
         const Bits& x_ns_nonred_in_y_closure) {
      for (auto i : index_range(x.cs.sk_rows)) {
        const auto& c_i = x.cs.sk_rows[i];
        // make sure c_i violates y
        if (y.relation_with(c_i).implies(Poly_Con_Rel::is_included()))
          continue;
        if (not x_ln_nonred_in_y.empty())
          // C1: c_i saturates a line non-redundant in y
          return true;
        // get generators saturated by c
        Bits sat_i = x.sat_g[i];
        sat_i.complement_until(x.sat_g.num_cols());
        if (sat_i.intersects(x_sk_np_nonred_in_y))
          // C1: c_i saturates a skel non-point non-redundant in y
          return true;
        if (c_i.is_nonstrict_inequality()) {
          if (sat_i.intersects(x_sk_p_nonred_in_y_closure))
            // C2: c_i saturates a skel point non-redundant in y closure
            return true;
          for (auto j : x_ns_nonred_in_y_closure) {
            const auto& ns_j = x.gs.ns_rows[j];
            if (ns_j.test(i))
              // C2: c_i belongs to the support (i.e., it saturates)
              // a non-skel point non-redundant in y closure
              return true;
          }
        }
      }
      return false;
    };

  auto check_C3_on_poly
    = [](const Poly_Impl& x, const Poly_Impl& y, const Poly_Impl& hull) {

      // helper: checks condition C3 on a single constraint;
      // returns true if union is NOT exact.
      auto check_C3 = [](const Con& c,
                         const Poly_Impl& y, const Poly_Impl& hull) {
        assert(c.is_strict_inequality());
        // create an equality from c
        auto eq = Con(c.linear_expr(), c.inhomo_term(), Con::EQUALITY);
        // interesect both y and hull with eq
        auto y_eq = y;
        y_eq.add_con(eq);
        auto hull_eq = hull;
        hull_eq.add_con(eq);
        // if y_eq does not contain hull_eq, hull is not exact
        return not y_eq.contains(hull_eq);
      };

      // first check skel constraints (no need to check singular ones)
      for (const auto& c : x.cs.sk_rows) {
        // make sure c is a strict ineq that violates y
        if (not c.is_strict_inequality()
            || y.relation_with(c).implies(Poly_Con_Rel::is_included()))
          continue;
        if (check_C3(c, y, hull))
          return true;
      }
      // then check non-skel constraints
      for (const auto& ns : x.cs.ns_rows) {
        const auto& c = detail::materialize(ns, x.cs.sk_rows);
        assert(c.is_strict_inequality());
        // make sure c violates y
        if (y.relation_with(c).implies(Poly_Con_Rel::is_included()))
          continue;
        if (check_C3(c, y, hull))
          return true;
      }
      return false;
    };


  // check x equalities wrt y (and vice versa)
  if (check_equalities(x, y,
                       x_ln_nonred_in_y, x_sk_np_nonred_in_y,
                       x_sk_p_nonred_in_y_closure, x_ns_nonred_in_y_closure))
    // hull not exact due to equality in x
    return false;
  if (check_equalities(y, x,
                       y_ln_nonred_in_x, y_sk_np_nonred_in_x,
                       y_sk_p_nonred_in_x_closure, y_ns_nonred_in_x_closure))
    // hull not exact due to equality in y
    return false;

  // check properties C1 and C2 on x skel ineqs wrt y (and vice versa)
  if (check_C1_C2_on_skel_ineqs(x, y,
                                x_ln_nonred_in_y, x_sk_np_nonred_in_y,
                                x_sk_p_nonred_in_y_closure,
                                x_ns_nonred_in_y_closure))
    // hull not exact due to skel ineq in x
    return false;
  if (check_C1_C2_on_skel_ineqs(y, x,
                                y_ln_nonred_in_x, y_sk_np_nonred_in_x,
                                y_sk_p_nonred_in_x_closure,
                                y_ns_nonred_in_x_closure))
    // hull not exact due to skel ineq in y
    return false;

  // compute the poly hull, needed to check C3.
  Poly_Impl hull(x);
  for (auto i : y_ln_nonred_in_x)
    hull.add_gen(y.gs.sing_rows[i]);
  for (auto i : y_sk_nonred_in_x)
    hull.add_gen(y.gs.sk_rows[i]);
  for (auto i : y_ns_nonred_in_x)
    hull.add_gen(y_gs_ns[i]);
  hull.minimize();

  // check for C3 on x strict ineqs wrt y and hull (and vice versa)
  if (check_C3_on_poly(x, y, hull))
    return false;
  if (check_C3_on_poly(y, x, hull))
    return false;

  // All checks passed: hull is exact.
  x = std::move(hull);
  return true;
}

} // namespace pplite
