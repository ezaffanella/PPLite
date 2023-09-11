/* PPLite: a lightweight library for convex polyhedra derived from PPL.
   Copyright (C) 2001-2010 Roberto Bagnara <bagnara@cs.unipr.it>
   Copyright (C) 2010-2018 BUGSENG srl (http://bugseng.com)
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

#include "pplite-config.h"
#include "utils.hh"
#include "Bits.hh"
#include "Con.hh"
#include "Gen.hh"
#include "Integer.hh"
#include "Linear_Expr.hh"
#include "Poly.hh"
#include "Poly_min.hh"
#include "Poly_widen.hh"
#include "Poly_Rel.hh"
#include "Sat.hh"
#include "Scalar_Prod.hh"

#include <algorithm>
#include <cassert>
#include <deque>
#include <functional>
#include <iterator>
#include <string>
#include <tuple>
#include <utility>

namespace pplite {

namespace h79_widen {

std::tuple<Index_Set, Index_Set>
select_risky_h79_unstable_cons(const Poly_Impl& x, const Poly_Impl& y) {
  assert(x.topol == y.topol && x.dim > 0 && x.dim == y.dim);
  assert(x.marked_min() && y.marked_min());
  const auto& y_sk = y.cs.sk_rows;
  // Detect the positivity constraint in y_sk (if present):
  // we widen the polyhedron, not the underlying polyhedral cone.
  auto pos_iter = std::find_if(y_sk.begin(), y_sk.end(),
                               std::mem_fn(&Con::is_tautological));
  const bool has_pos = (pos_iter != y_sk.end());

  // *Copy* y.sat_g rows.
  auto s_rows = y.sat_g.impl().rows;
  if (has_pos) {
    // Remove it
    auto pos_idx = std::distance(y_sk.begin(), pos_iter);
    using std::swap;
    swap(s_rows[pos_idx], s_rows.back());
    s_rows.pop_back();
  }
  // Sort sat rows
  std::sort(s_rows.begin(), s_rows.end());

  Index_Set sk_unstable;
  Index_Set ns_unstable;
  // Process skel constraints.
  Bits buffer;
  for (auto i : bwd_index_range(x.cs.sk_rows)) {
    const auto& c_i = x.cs.sk_rows[i];
    if (c_i.is_tautological())
      // c_i is stable: positivity.
      continue;
    buffer.clear();
    for (auto j : bwd_index_range(y.gs.sk_rows)) {
      const auto& g_j = y.gs.sk_rows[j];
      const int sp_sgn = sp::sign(c_i, g_j);
      // We are assuming that `y <= x'.
      assert(sp_sgn >= 0);
      if (sp_sgn > 0)
        buffer.set(j);
    }
    if (buffer.empty())
      // c_i is stable: same saturators of an equality of y.
      continue;
    if (std::binary_search(s_rows.begin(), s_rows.end(), buffer))
      // c_i is stable: same saturators of a skel ineq of y.
      continue;
    // c_i is unstable.
    sk_unstable.set(i);
  }
  Index_Set sk_stable = sk_unstable;
  sk_stable.complement_until(num_rows(x.cs.sk_rows));
  // Process non-skel constraints.
  for (auto i : bwd_index_range(x.cs.ns_rows)) {
    const auto& ns_i = x.cs.ns_rows[i];
    // If the whole support is stable, ns_i is stable too.
    if (subset_eq(ns_i, sk_stable))
      continue;
    ns_unstable.set(i);
  }
  return std::make_tuple(sk_unstable, ns_unstable);
}

void
risky_h79_widen(Poly_Impl& x, const Poly_Impl& y) {
  assert(x.topol == y.topol);
  assert(x.dim > 0 && x.dim == y.dim);
  assert(x.marked_min() && y.marked_min());

  Index_Set sk_unstable;
  Index_Set ns_unstable;
  std::tie(sk_unstable, ns_unstable) = select_risky_h79_unstable_cons(x, y);
  // If no skel constraint was unstable, the result is `x'.
  // (Note: equalities are always stable, since x contains y;
  // also, if all skel ineqs are stable, then all non-skel are stable too.)
  // The test allows for preserving minimal form.
  if (sk_unstable.empty())
    return;
  // We found some unstable constraints: drop them.
  Poly_Impl::Sys<Cons> h79_cs = std::move(x.cs);
  erase_using_sorted_indices(h79_cs.sk_rows, sk_unstable);
  erase_using_sorted_indices(h79_cs.ns_rows, ns_unstable);
  // Update ns indices.
  for (auto& ns : h79_cs.ns_rows)
    ns.remove_all(sk_unstable);
  x.set_universe();
  if (!h79_cs.empty()) {
    x.cs_pending = std::move(h79_cs);
    x.set_status(Poly_Impl::Status::PENDING);
  }
  assert(x.check_inv());
}

std::tuple<Index_Set, Index_Set>
select_safe_h79_unstable_cons(const Poly_Impl& x, const Poly_Impl& y) {
  assert(x.topol == y.topol && x.dim > 0 && x.dim == y.dim);
  assert(x.marked_min() && y.marked_min());

  Index_Set sk_unstable;
  Index_Set ns_unstable;
  // Process y's skel constraints.
  const auto& y_sk = y.cs.sk_rows;
  for (auto i : bwd_index_range(y_sk)) {
    if (not sp::satisfied_by_all(x.topol, y_sk[i], x.gs))
      sk_unstable.set(i);
  }
  // Process non-skel constraints.
  Index_Set sk_stable = sk_unstable;
  sk_stable.complement_until(num_rows(y_sk));
  const auto& y_ns = y.cs.ns_rows;
  for (auto i : bwd_index_range(y_ns)) {
    // If the whole support is stable, ns_i is stable too.
    if (subset_eq(y_ns[i], sk_stable))
      continue;
    ns_unstable.set(i);
  }
  return std::make_tuple(sk_unstable, ns_unstable);
}

void
safe_h79_widen(Poly_Impl& x, const Poly_Impl& y) {
  assert(x.topol == y.topol);
  assert(x.dim > 0 && x.dim == y.dim);
  assert(x.marked_min() && y.marked_min());

  // Current implementation only works for closed polyhedra.
  assert(x.topol == Topol::CLOSED);

  Index_Set sk_unstable;
  Index_Set ns_unstable;
  std::tie(sk_unstable, ns_unstable) = select_safe_h79_unstable_cons(x, y);
  if (sk_unstable.empty()) {
    x = y;
    return;
  }
  // We found some unstable constraints: drop them.
  Poly_Impl::Sys<Cons> h79_cs = y.cs;
  erase_using_sorted_indices(h79_cs.sk_rows, sk_unstable);
  erase_using_sorted_indices(h79_cs.ns_rows, ns_unstable);

  // Update ns indices.
  for (auto& ns : h79_cs.ns_rows)
    ns.remove_all(sk_unstable);
  x.set_universe();
  if (!h79_cs.empty()) {
    x.cs_pending = std::move(h79_cs);
    x.set_status(Poly_Impl::Status::PENDING);
  }
  assert(x.check_inv());
}

void
h79_widen(Poly_Impl& x, const Poly_Impl& y, Widen_Spec w_spec) {
  if (w_spec == Widen_Spec::SAFE) {
    if (x.topol == Topol::CLOSED)
      safe_h79_widen(x, y);
    else {
      // For NNC, use trivial lifting of risky widening.
      x.poly_hull_assign(y);
      x.minimize();
      risky_h79_widen(x, y);
    }
  } else
    risky_h79_widen(x, y);
}

void
boxed_h79_widen(Poly_Impl& x, const Poly_Impl& y, Widen_Spec w_spec) {
  assert(x.topol == y.topol);
  assert(x.dim > 0 && x.dim == y.dim);
  assert(x.marked_min() && y.marked_min());

  auto x_bbox = x.get_bounding_box();
  const auto& y_bbox = y.get_bounding_box();
  x_bbox.widening_assign(y_bbox);

  h79_widen(x, y, w_spec);

  Cons boxed_cs;
  for (auto i : dim_range(x_bbox)) {
    Var var(i);
    const auto& itv = x_bbox.get_bounds(var);
    if (itv.is_singleton()) {
      assert(x.relation_with(get_eq_con(var, itv))
             .implies(Poly_Con_Rel::is_included()));
      continue;
    }
    if (itv.has_lb())
      boxed_cs.push_back(get_lb_con(var, itv));
    if (itv.has_ub())
      boxed_cs.push_back(get_ub_con(var, itv));
  }
  x.add_cons(boxed_cs);
  assert(x.check_inv());
}

} // namespace h79_widen

namespace bhrz03_widen {

Cert::Cert(const Poly_Impl& ph)
  : num_sk_points(0), num_rays_null_coord(ph.dim, 0) {
  // Need non-redundant descriptions.
  ph.minimize();
  // Precondition: `ph' is not empty.
  assert(ph.marked_min());

  affine_dim = ph.affine_dim();
  lin_space_dim = ph.num_lines();

#ifndef NDEBUG
  equalities = ph.num_equals();
#endif
  // Count skel constraints, except for the positivity.
  // While at it, check if the positivity constraint is strict.
  const auto& sk_cons = ph.cs.sk_rows;
  num_sk_cons = num_rows(sk_cons);
  const auto pos_iter = std::find_if(sk_cons.cbegin(), sk_cons.cend(),
                                     std::mem_fn(&Con::is_tautological));
  if (pos_iter != sk_cons.cend()) {
    --num_sk_cons;
    assert(pos_iter->is_strict_inequality());
  }

  if (ph.topol == Topol::NNC) {
    auto sk_strict = std::count_if(sk_cons.cbegin(), sk_cons.cend(),
                                   std::mem_fn(&Con::is_strict_inequality));
    if (sk_strict > 0)
      support_cards[1] = sk_strict;
    for (const auto& ns : ph.cs.ns_rows)
      ++support_cards[ns.size()];
    // Check for the efc (and erase it, if found).
    Index_Set rfc = detail::find_efc(ph.cs.ns_rows,
                                     ph.cs.sk_rows, ph.gs.sk_rows,
                                     ph.sat_g, ph.sat_c);
    const auto rfc_size = rfc.size();
    if (rfc_size > 0) {
      assert(support_cards[rfc_size] > 0);
      assert(detail::is_empty_face_cutter(rfc, ph.sat_g,
                                          ph.gs.sk_rows));
      --support_cards[rfc_size];
      if (support_cards[rfc_size] == 0)
        support_cards.erase(rfc_size);
    }
  }

  assert(num_sk_points == 0);
  for (const auto& g : ph.gs.sk_rows) {
    if (g.is_ray()) {
      // For each z such that 0 <= z < ph.dim,
      // `num_rays_null_coord[z]' will be the number of rays
      // having exactly `z' coordinates equal to 0.
      auto zeroes = g.linear_expr().num_zeroes(0, g.space_dim());
      zeroes += (ph.dim - g.space_dim());
      ++num_rays_null_coord[zeroes];
    } else {
      assert(g.is_point() || g.is_closure_point());
      ++num_sk_points;
    }
  }
  assert(check_inv());
}

Cert::Cert(const F_Poly& ph) {
  ph.minimize();
  assert(not ph.is_empty());
  const auto sd = ph.space_dim();
  num_rays_null_coord.resize(sd, 0);

  // Scan interval factors.
  for (const auto& itv : ph.impl().proper_itvs()) {
    assert(not itv.is_empty());
    if (itv.is_universe()) {
      ++affine_dim;
      ++lin_space_dim;
    } else if (itv.is_bounded()) {
      if (not itv.is_singleton()) {
        ++affine_dim;
        num_sk_cons += 2;
        num_sk_points *= 2;
      } else {
        // Singleton
#ifndef NDEBUG
        ++equalities;
#endif
      }
    } else {
      // Not bounded, not universe.
      ++affine_dim;
      ++num_sk_cons;
      ++num_rays_null_coord[sd-1];
    }
  }

  // Scan polyhedra factors
  const auto& fs = ph.impl().factors;
  const auto nf = num_rows(fs);
  for (auto i : range(nf)) {
    const auto& f_i = fs[i];
    Cert cert_i(f_i.impl());

    affine_dim += cert_i.affine_dim;
    lin_space_dim += cert_i.lin_space_dim;
    num_sk_cons += cert_i.num_sk_cons;
#ifndef NDEBUG
    equalities += cert_i.equalities;
#endif
    for (const auto& p : cert_i.support_cards)
      support_cards[p.first] += p.second;
    num_sk_points *= cert_i.num_sk_points;
    const auto sd_i = f_i.space_dim();
    const auto delta = sd - sd_i;
    for (auto d : range(sd_i))
      num_rays_null_coord[delta + d] += cert_i.num_rays_null_coord[d];
  }
  assert(check_inv());
#ifndef NDEBUG
  Poly ph_poly = ph.to_poly();
  Cert cert_poly(ph_poly.impl());
  assert(compare(cert_poly) == 0);
#endif
}

int
Cert::compare(const Cert& y) const {
  const auto& x = *this;
  assert(x.check_inv() && y.check_inv());
  if (x.affine_dim != y.affine_dim)
    return (x.affine_dim < y.affine_dim) ? 1 : -1;
  if (x.lin_space_dim != y.lin_space_dim)
    return (x.lin_space_dim < y.lin_space_dim) ? 1 : -1;
  if (x.num_sk_cons != y.num_sk_cons)
    return (x.num_sk_cons > y.num_sk_cons) ? 1 : -1;
  { // Comparing support-size cardinalities.
    auto x_i = x.support_cards.cbegin(), x_e = x.support_cards.cend();
    auto y_i = y.support_cards.cbegin(), y_e = y.support_cards.cend();
    for ( ; x_i != x_e && y_i != y_e; ++x_i, ++y_i) {
      auto x_size = x_i->first;
      auto y_size = y_i->first;
      if (x_size != y_size)
        // The smaller of the two support sizes can not be found
        // in the other polyhedron; when the skeleton is invariant,
        // a smaller support size corresponds to a bigger affine
        // dimension for the cut face.
        return (x_size < y_size) ? 1 : -1;
      auto x_card = x_i->second;
      auto y_card = y_i->second;
      if (x_card != y_card)
        // Different cardinalities for the given support size;
        // when the skeleton is invariant, this means the polyhedron
        // have a different number of cut faces of a specific affine dimension.
        return (x_card > y_card) ? 1 : -1;
    }
    if (x_i != x_e || y_i != y_e)
      // Ditto: some support sizes can not be found in the other polyhedron.
      return (x_i != x_e) ? 1 : -1;
  }
  if (x.num_sk_points != y.num_sk_points)
    return (x.num_sk_points > y.num_sk_points) ? 1 : -1;
  { // Comparing rays null coords.
    const dim_type dims = x.num_rays_null_coord.size();
    assert(x.num_rays_null_coord.size() == y.num_rays_null_coord.size());
    // Note: iterating upwards, because we have to check first
    // the number of rays having more NON-zero coordinates.
    for (auto i : range(dims)) {
      if (x.num_rays_null_coord[i] != y.num_rays_null_coord[i])
        return (x.num_rays_null_coord[i] > y.num_rays_null_coord[i]) ? 1 : -1;
    }
  }
  // All components are equal.
  return 0;
}

int
Cert::compare(const Poly_Impl& ph) const {
  Cert ph_cert(ph);
  return compare(ph_cert);
}

bool
Cert::check_inv() const {
  std::string reason;
#ifdef NDEBUG
  auto maybe_dump = []() {};
#else // In debugging mode, be noisy.
  auto maybe_dump = [&reason]() {
    std::cerr << reason << std::endl;
  };
#endif

  // The dimension of the vector space.
  const dim_type space_dim = num_rays_null_coord.size();
  if (affine_dim > space_dim) {
    reason = "bhrz03 cert: affine dim greater than space dim";
    maybe_dump();
    return false;
  }
  if (lin_space_dim > affine_dim) {
    reason = "bhrz03 cert: lineality space dim greater than affine dim";
    maybe_dump();
    return false;
  }
#ifndef NDEBUG
  if (num_sk_cons + equalities < space_dim - affine_dim) {
    reason = "bhrz03 cert: too few constraints";
    maybe_dump();
    return false;
  }
#endif
  if (num_sk_points == 0) {
    reason = "bhrz03 cert: no (closure) points";
    maybe_dump();
    return false;
  }
  if (lin_space_dim == space_dim) {
    if (num_sk_cons > 0) {
      reason = "bhrz03 cert: universe poly has non-redundant constraints";
      maybe_dump();
      return false;
    }
    if (num_sk_points != 1) {
      reason = "bhrz03 cert: universe poly has more than one point";
      maybe_dump();
      return false;
    }
  }

  // All tests passed.
  return true;
}

bool
combining_cons(Poly_Impl& x, const Poly_Impl& y, const Cert& y_cert,
               const Poly_Impl& h79, const Index_Set& sk_unstable) {
  // We will choose from `sk_unstable' many subsets of constraints.
  // For each subset, an average constraint will be added to `new_cs'.
  // There is no point in applying this technique when `sk_unstable'
  // has only one element.
  if (sk_unstable.size() <= 1)
    return false;

  Cons new_cs;
  // Consider the skel (closure) points that belong to
  // the generator systems of both x and y.
  const auto& h79_cs_sk = h79.cs.sk_rows;
  for (const auto& g : y.gs.sk_rows) {
    if (g.is_ray())
      continue;
    // If in the skel of `h79' there is already an inequality
    // saturating this (closure) point, no need to produce a constraint.
    if (any_of(h79_cs_sk, [&g](const Con& c) {
                           return sp::sign(c, g) == 0;
                          }))
      continue;
    // Consider all the constraints in `x_cs_unstable_sk'
    // that are saturated by the point `g'.
    bool first_con = true;
    Linear_Expr expr;
    Integer inhomo;
    Con::Type type = Con::NONSTRICT_INEQUALITY;
    for (dim_type i : sk_unstable) {
      const auto& c = x.cs.sk_rows[i];
      if (sp::sign(c, g) != 0)
        continue;
      if (first_con) {
        first_con = false;
        expr = c.linear_expr();
        inhomo = c.inhomo_term();
        if (c.is_strict_inequality())
          type = Con::STRICT_INEQUALITY;
        continue;
      }
      // Not the first: combine constraints.
      expr += c.linear_expr();
      inhomo += c.inhomo_term();
      if (c.is_strict_inequality())
        type = Con::STRICT_INEQUALITY;
    }
    if (expr.is_zero())
      continue;
    new_cs.emplace_back(expr, inhomo, type);
  }

  // Check for improvement wrt `h79'.
  if (none_of(new_cs, [&h79](const Con& c) {
                        return h79.relation_with(c)
                          == Poly_Con_Rel::strictly_intersects();
                      }))
    return false;

  Poly_Impl result = h79;
  result.add_cons(std::move(new_cs));
  // Force minimization.
  result.minimize();
  if (!y_cert.is_stabilizing(result))
    return false;
  // Heuristics was successful.
  x = std::move(result);
  return true;
}

bool
evolving_points(Poly_Impl& x, const Poly_Impl& y, const Cert& y_cert,
                const Poly_Impl& h79) {
  // For each skel (closure) point in `x' that is not in `y',
  // this technique tries to identify a set of rays that:
  //  - are included in polyhedron `h79';
  //  - when added to `y' will subsume the point.
  Gens cand_rays;
  for (const auto& x_g : x.gs.sk_rows) {
    if (x_g.is_ray())
      continue;
    if (y.relation_with(x_g) == Poly_Gen_Rel::subsumes())
      continue;
    // For each skel (closure) point `y_g' in `y',
    // build the candidate ray `x_g - y_g'.
    for (const auto& y_g : y.gs.sk_rows) {
      if (y_g.is_ray())
        continue;
      assert(compare(x_g, y_g) != 0);
      // This is a variant of Gen::linear_combine,
      // but we zero the divisor rather than a space dimension.
      auto x_e = x_g.linear_expr();
      auto x_d = x_g.divisor();
      Linear_Expr::combine(x_e, x_d,
                           y_g.linear_expr(), y_g.divisor(),
                           x_g.divisor(), y_g.divisor());
      assert(x_d == 0 && !x_e.is_zero());
      cand_rays.push_back(ray(x_e));
    }
  }
  if (cand_rays.empty())
    return false;

  Poly_Impl result = x;
  result.add_gens(std::move(cand_rays));
  result.intersection_assign(h79);
  result.minimize();
  if (!y_cert.is_stabilizing(result) || result.contains(h79))
    return false;
  // Heuristics was successful.
  x = std::move(result);
  return true;
}

Gen
evolve_ray(const Gen& x_ray, const Gen& y_ray) {
  const auto& x_expr = x_ray.linear_expr();
  const auto& y_expr = y_ray.linear_expr();
  Linear_Expr ray_expr = x_expr;

  const auto x_sd = x_ray.space_dim();
  std::deque<bool> considered(x_sd, false);
  Integer tmp;
  for (auto i : range(x_sd)) {
    if (x_expr[i].is_zero() || considered[i])
      continue;
    for (auto j : range(i+1, x_sd)) {
      if (x_expr[j].is_zero() || considered[j])
        continue;
      tmp = x_expr[i] * y_expr.get(j);
      sub_mul_assign(tmp, x_expr[j], y_expr.get(i));
      const int clockwise = sgn(tmp);
      const int first_or_third_quadrant = sgn(x_expr[i]) * sgn(x_expr[j]);
      switch (clockwise * first_or_third_quadrant) {
      case -1:
        ray_expr.set(i, Integer::zero());
        considered[i] = true;
        break;
      case 1:
        ray_expr.set(j, Integer::zero());
        considered[j] = true;
        break;
      default:
        break;
      }
    }
  }
  return ray(ray_expr);
}

bool
evolving_rays(Poly_Impl& x, const Poly_Impl& y, const Cert& y_cert,
              const Poly_Impl& h79) {
  Gens cand_rays;
  for (const auto& x_g : x.gs.sk_rows) {
    // Check that it is a ray not belonging to `y'.
    if (!x_g.is_ray() || y.relation_with(x_g) == Poly_Gen_Rel::subsumes())
      continue;
    for (const auto& y_g : y.gs.sk_rows) {
      if (y_g.is_ray())
        cand_rays.push_back(evolve_ray(x_g, y_g));
    }
  }
  if (cand_rays.empty())
    return false;

  Poly_Impl result = x;
  result.add_gens(std::move(cand_rays));
  result.intersection_assign(h79);
  result.minimize();
  if (!y_cert.is_stabilizing(result) || result.contains(h79))
    return false;
  // Heuristics was successful.
  x = std::move(result);
  return true;
}

void
risky_bhrz03_widen(Poly_Impl& x, const Poly_Impl& y) {
  assert(x.topol == y.topol);
  assert(x.dim > 0 && x.dim == y.dim);
  assert(x.marked_min() && y.marked_min());

  // Compute certificate info for `y'.
  const Cert y_cert(y);

  // If the iteration is stabilizing, the result is `x'.
  // Also check if the two polyhedra are the same
  // (exploiting the knowledge that `y <= x').
  if (y_cert.is_stabilizing(x) || y.contains(x))
    return;

  // Here the iteration is not immediately stabilizing.
  Index_Set sk_unstable;
  Index_Set ns_unstable;
  std::tie(sk_unstable, ns_unstable)
    = h79_widen::select_risky_h79_unstable_cons(x, y);
  // There has to be some unstable skel constraint, since otherwise
  // the iteration should have been immediately stabilizing.
  assert(!sk_unstable.empty());

  // Compute the h79 widening.
  // NOTE: we *copy* constraints, since we can not modify `x'.
  Poly_Impl::Sys<Cons> h79_cs;
  h79_cs.sing_rows = x.cs.sing_rows;
  for (auto i : index_range(x.cs.sk_rows)) {
    if (sk_unstable.test(i))
      continue;
    h79_cs.sk_rows.push_back(x.cs.sk_rows[i]);
  }
  for (auto i : index_range(x.cs.ns_rows)) {
    if (ns_unstable.test(i))
      continue;
    auto ns = x.cs.ns_rows[i];
    ns.remove_all(sk_unstable);
    h79_cs.ns_rows.push_back(std::move(ns));
  }

  Poly_Impl h79(x.dim, Spec_Elem::UNIVERSE, x.topol);
  if (!h79_cs.empty()) {
    h79.cs_pending = std::move(h79_cs);
    h79.set_status(Poly_Impl::Status::PENDING);
  }
  // Force minimization.
  h79.minimize();

  assert(x.topol == y.topol && x.topol == h79.topol);
  assert(x.dim > 0 && x.dim == y.dim && x.dim == h79.dim);
  assert(x.marked_min() && y.marked_min() && h79.marked_min());

  // None of the following widening heuristics is intrusive:
  // they will modify `x' only when returning successfully.
  if (combining_cons(x, y, y_cert, h79, sk_unstable))
    return;
  if (evolving_points(x, y, y_cert, h79))
    return;
  if (evolving_rays(x, y, y_cert, h79))
    return;

  // No previous technique was successful: fall back to the h79 widening.
  x = std::move(h79);
  assert(x.check_inv());
  // The h79 widening is always stabilizing.
  assert(y_cert.is_stabilizing(x));
}

void
bhrz03_widen(Poly_Impl& x, const Poly_Impl& y, Widen_Spec w_spec) {
  assert(x.topol == y.topol);
  assert(x.dim > 0 && x.dim == y.dim);
  assert(x.marked_min() && y.marked_min());

  if (w_spec == Widen_Spec::SAFE) {
    // Trivial lifting of risky spec.
    x.poly_hull_assign(y);
    x.minimize();
  }
  risky_bhrz03_widen(x, y);
}

} // namespace bhrz03_widen

} // namespace pplite


