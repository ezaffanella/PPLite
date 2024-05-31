/* PPLite: a lightweight library for convex polyhedra derived from PPL.
   Copyright (C) 2001-2010 Roberto Bagnara <bagnara@cs.unipr.it>
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

#include "pplite-config.h"
#include "utils.hh"
#include "Poly.hh"
#include "Bits.hh"
#include "Integer.hh"
#include "Scalar_Prod.hh"
#include "Poly_templ.hh"
#include "Poly_min.hh"
#include "Poly_widen.hh"
#include "Rational.hh"

#include <algorithm>
#include <functional>
#include <numeric> // for std::iota
#include <string>
#include <tuple>
#include <unordered_map>

namespace pplite {

namespace detail {

template <typename T1, typename T2>
struct Pair_Hasher {
  size_t operator()(const std::pair<T1, T2>& p) const {
    size_t res = std::hash<T1>{}(p.first);
    hash_combine(res, std::hash<T2>{}(p.second));
    return res;
  }
};

template <typename Row>
size_t simple_hash(const Row& row) {
  const auto sd = row.space_dim();
  const auto t = static_cast<size_t>(row.type());
  const auto& ex = row.impl().expr;
  const auto& inhomo = row.impl().inhomo;
  size_t res = hash_size(t);
  hash_combine(res, inhomo.hash());
  for (auto i : range(sd)) {
    if (ex[i] != 0)
      hash_combine(res, (1+i) * ex[i].hash());
  }
  return res;
}

template <typename Row>
size_t simple_hash(const std::vector<Row>& rows) {
  // We partially (and indirectly) sort the container:
  // 10 elements at most should be enough.
  const size_t max_sorted = 10;
  const size_t sz = rows.size();
  const size_t num_sorted = std::min(sz, max_sorted);

  auto indirect_cmp = [&rows](dim_type i, dim_type j) {
    return compare(rows[i], rows[j]) < 0;
  };

  Dims idxs(sz);
  std::iota(idxs.begin(), idxs.end(), 0);
  std::partial_sort(idxs.begin(), idxs.begin() + num_sorted,
                    idxs.end(), indirect_cmp);
  // Hash the size and the first num_sorted rows.
  size_t res = hash_size(sz);
  for (size_t i = 0; i < num_sorted; ++i)
    hash_combine(res, simple_hash(rows[idxs[i]]));
  // Also hash last row, if not already done.
  if (num_sorted < sz) {
    auto iter = std::max_element(idxs.begin() + num_sorted, idxs.end(),
                                 indirect_cmp);
    hash_combine(res, simple_hash(rows[*iter]));
  }
  return res;
}

} // namespace detail

Poly_Impl::Cons_Proxy
Poly_Impl::normalized_cons() const {
  detail::normalize_for_hash(cs);
  return cons();
}

size_t
Poly_Impl::hash() const {
  minimize();
  using detail::normalize_for_hash;
  using detail::simple_hash;
  size_t res = hash_size(dim);
  if (is_empty())
    return res;
  // Normalize and hash constraint system.
  normalize_for_hash(cs);
  hash_combine(res, simple_hash(cs.sing_rows));
  hash_combine(res, simple_hash(cs.sk_rows));
  // For ns, hash only the size.
  hash_combine(res, hash_size(cs.ns_rows.size()));
  // Note: avoid normalizing generator system,
  // as it often has many elements; just hash the sizes.
  hash_combine(res, hash_size(gs.sing_rows.size()));
  hash_combine(res, hash_size(gs.sk_rows.size()));
  hash_combine(res, hash_size(gs.ns_rows.size()));
  return res;
}

namespace detail {

bool
check_sat(const Sat& sat_g, const Cons& cons, const Gens& gens) {
  for (auto i : bwd_index_range(cons)) {
    const auto& sat_i = sat_g[i];
    for (auto j : bwd_index_range(gens)) {
      const bool saturated = !sat_i[j];
      const int s = sp::sign(cons[i], gens[j]);
      if (s < 0) return false;
      if (saturated && s > 0) return false;
      if (!saturated && s == 0) return false;
    }
  }
  return true;
}

} // namespace detail

void
Poly_Impl::clear_all() {
  // Clears everything but topol, dim and status.
  cs.clear();
  gs.clear();
  sat_c.clear();
  sat_g.clear();
  cs_pending.clear();
  gs_pending.clear();
}

void
Poly_Impl::set_empty() {
  clear_all();
  // Add `1 == 0'.
  cs.sing_rows.push_back(Con::zero_dim_false());
  set_status(Status::EMPTY);
  assert(check_inv());
}

void
Poly_Impl::set_universe() {
  clear_all();
  detail::set_universe(dim, cs, gs, sat_c);
  sat_g = sat_c.transpose();
  set_status(Status::MINIMIZED);
  // check_inv() would trigger an infinite recursion.
  assert(check_inv_without_minimize());
}

void
Poly_Impl::set_topology(Topol t) {
  // Topological closure is assumed to hold.
  assert(topol == Topol::CLOSED || t == Topol::NNC
         || is_topologically_closed());
  topol = t;
}

void
Poly_Impl::reinit_with_gens(dim_type d, Topol t, Sys<Gens>&& gens) {
  assert(&gens != &gs && &gens != &gs_pending);
  dim = d;
  topol = t;
  clear_all();
  gs_pending = std::move(gens);
  detail::init_dd(dim, gs, gs_pending, cs, sat_g);
  sat_c = sat_g.transpose();
  set_status(gs_pending.empty() ? Status::MINIMIZED : Status::PENDING);
  // check_inv() would trigger an infinite recursion.
  assert(check_inv_without_minimize());
}

Poly_Impl::Poly_Impl(dim_type d, Spec_Elem s, Topol t)
  : topol(t), dim(d) {
  if (s == Spec_Elem::EMPTY)
    set_empty();
  else
    set_universe();
}

bool
Poly_Impl::is_empty() const {
  if (has_valid_gens()) {
    // (NS_g = \empty && SK_g has no point) --> marked_empty()
    assert(!(gs.ns_rows.empty() && !has_point(gs.sk_rows))
           || marked_empty());
    return marked_empty();
  }
  minimize();
  return marked_empty();
}

bool
Poly_Impl::is_universe() const {
  if (marked_empty())
    return false;
  if (space_dim() == 0)
    return true;
  if (has_valid_cons()) {
    auto& pend_sing = cs_pending.sing_rows;
    auto& pend_sk = cs_pending.sk_rows;
    auto& pend_ns = cs_pending.ns_rows;
    const auto is_tauto = std::mem_fn(&Con::is_tautological);
    return cs.sing_rows.empty() && cs.ns_rows.empty()
      && all_of(cs.sk_rows, is_tauto)
      && pend_sing.empty() && pend_ns.empty()
      && all_of(pend_sk, is_tauto);
  }
  minimize();
  return !marked_empty() && (num_lines() == space_dim());
}

bool
Poly_Impl::is_bounded() const {
  if (marked_empty() || space_dim() == 0)
    return true;
  if (has_valid_gens())
    return gs.sing_rows.empty()
      && !has_line_or_ray(gs.sk_rows)
      && gs_pending.sing_rows.empty()
      && !has_line_or_ray(gs_pending.sk_rows);
  minimize();
  return gs.sing_rows.empty() && !has_line_or_ray(gs.sk_rows);
}

bool
Poly_Impl::is_topologically_closed() const {
  if (is_necessarily_closed())
    return true;
  if (marked_empty() || space_dim() == 0)
    return true;
  minimize();
  return gs.ns_rows.empty() && !has_closure_point(gs.sk_rows);
}

dim_type
Poly_Impl::affine_dim() const {
  minimize();
  return marked_empty() ? 0 : space_dim() - cs.sing_rows.size();
}

dim_type
Poly_Impl::num_min_cons() const {
  minimize();
  auto num_cons = cs.sing_rows.size() + cs.sk_rows.size() + cs.ns_rows.size();
  bool has_pos = any_of(cs.sk_rows, std::mem_fn(&Con::is_tautological));
  if (has_pos)
    return num_cons - 1;
  bool has_efc = detail::has_efc(cs.ns_rows, sat_g, gs.sk_rows);
  if (has_efc)
    return num_cons - 1;
  // No positivity nor efc.
  return num_cons;
}

Gens_Info
Poly_Impl::gens_info() const {
  minimize();
  if (marked_empty())
    return { 0, 0, 0, 0, 0 };

  auto ln = num_rows(gs.sing_rows);
  auto r = 0, cp = 0, skp = 0;
  auto ns = num_rows(gs.ns_rows);

  for (const auto& g : gs.sk_rows) {
    switch (g.type()) {
    case Gen::POINT:
      ++skp;
      break;
    case Gen::RAY:
      ++r;
      break;
    case Gen::CLOSURE_POINT:
      ++cp;
      break;
    case Gen::LINE:
    default:
      PPLITE_UNREACH;
      break;
    }
  }
  return { ln, r, cp, skp, ns };
}

dim_type
Poly_Impl::num_min_gens() const {
  minimize();
  return num_rows(gs.sing_rows) + num_rows(gs.sk_rows) + num_rows(gs.ns_rows);
}

namespace detail {

Poly_Con_Rel
relation_with(const Poly_Impl::Gens_Proxy& gs, const Con& c) {
  Poly_Con_Rel res = Poly_Con_Rel::saturates();

  switch (c.type()) {

  case Con::EQUALITY:
    {
      res = res && Poly_Con_Rel::is_included();
      // The following variable will hold the scalar product sign
      // of either the first point or the first non-saturating ray.
      // It starts equal to 2, meaning that we haven't found any yet.
      int first_point_or_nonsat_ray_sign = 2;
      for (const auto& g : gs) {
        const int sp_sign = sp::sign(c, g);
        if (sp_sign == 0) {
          if (g.is_point()) {
            if (first_point_or_nonsat_ray_sign == 2)
              first_point_or_nonsat_ray_sign = 0;
            else if (first_point_or_nonsat_ray_sign != 0)
              return Poly_Con_Rel::strictly_intersects();
          }
        } else {
          assert(sp_sign != 0);
          if (g.is_line())
            return Poly_Con_Rel::strictly_intersects();
          // Non-saturating rays and closure points are treated as points.
          if (sp_sign == first_point_or_nonsat_ray_sign)
            continue;
          if (first_point_or_nonsat_ray_sign == 2) {
            first_point_or_nonsat_ray_sign = sp_sign;
            res = Poly_Con_Rel::is_disjoint();
          } else
            return Poly_Con_Rel::strictly_intersects();
        }
      }
    }
    break;

  case Con::NONSTRICT_INEQUALITY:
    {
      res = res && Poly_Con_Rel::is_included();
      bool no_point_or_nonsat_ray = true;
      for (const auto& g : gs) {
        const int sp_sign = sp::sign(c, g);
        if (sp_sign == 0) {
          if (g.is_point()) {
            if (no_point_or_nonsat_ray)
              no_point_or_nonsat_ray = false;
            else if (res == Poly_Con_Rel::is_disjoint())
              return Poly_Con_Rel::strictly_intersects();
          }
        } else {
          assert(sp_sign != 0);
          switch (g.type()) {
          case Gen::LINE:
            return Poly_Con_Rel::strictly_intersects();
          case Gen::RAY:
            if (no_point_or_nonsat_ray) {
              no_point_or_nonsat_ray = false;
              res = (sp_sign > 0)
                ? Poly_Con_Rel::is_included()
                : Poly_Con_Rel::is_disjoint();
            } else {
              if ((sp_sign > 0 && res == Poly_Con_Rel::is_disjoint())
                  ||
                  (sp_sign < 0 && res.implies(Poly_Con_Rel::is_included())))
                return Poly_Con_Rel::strictly_intersects();
              if (sp_sign > 0)
                res = Poly_Con_Rel::is_included();
            }
            break;
          case Gen::POINT:
          case Gen::CLOSURE_POINT:
            // A non-saturating closure point is treated as a point.
            if (no_point_or_nonsat_ray) {
              no_point_or_nonsat_ray = false;
              if (sp_sign > 0)
                res = Poly_Con_Rel::is_included();
              else if (sp_sign < 0)
                res = Poly_Con_Rel::is_disjoint();
            } else {
              if ((sp_sign > 0 && res == Poly_Con_Rel::is_disjoint())
                  ||
                  (sp_sign < 0 && res.implies(Poly_Con_Rel::is_included())))
                return Poly_Con_Rel::strictly_intersects();
              if (sp_sign > 0)
                res = Poly_Con_Rel::is_included();
            }
            break;
          } // switch (g.type())
        }
      }
    }
    break;

  case Con::STRICT_INEQUALITY:
    {
      res = res && Poly_Con_Rel::is_disjoint();
      bool no_point_or_nonsat_ray = true;
      for (const auto& g : gs) {
        const int sp_sign = sp::sign(c, g);
        if (sp_sign == 0) {
          if (g.is_point()) {
            if (no_point_or_nonsat_ray)
              no_point_or_nonsat_ray = false;
            else if (res == Poly_Con_Rel::is_included())
              return Poly_Con_Rel::strictly_intersects();
          }
        } else {
          assert(sp_sign != 0);
          switch (g.type()) {
          case Gen::LINE:
            return Poly_Con_Rel::strictly_intersects();
          case Gen::RAY:
            if (no_point_or_nonsat_ray) {
              no_point_or_nonsat_ray = false;
              res = (sp_sign > 0)
                ? Poly_Con_Rel::is_included()
                : Poly_Con_Rel::is_disjoint();
            } else {
              if ((sp_sign > 0 && res.implies(Poly_Con_Rel::is_disjoint()))
                  ||
                  (sp_sign <= 0 && res == Poly_Con_Rel::is_included()))
                return Poly_Con_Rel::strictly_intersects();
              if (sp_sign < 0)
                res = Poly_Con_Rel::is_disjoint();
            }
            break;
          case Gen::POINT:
          case Gen::CLOSURE_POINT:
            if (no_point_or_nonsat_ray) {
              no_point_or_nonsat_ray = false;
              if (sp_sign > 0)
                res = Poly_Con_Rel::is_included();
              else if (sp_sign < 0)
                res = Poly_Con_Rel::is_disjoint();
            } else {
              if ((sp_sign > 0 && res.implies(Poly_Con_Rel::is_disjoint()))
                  ||
                  (sp_sign <= 0 && res == Poly_Con_Rel::is_included()))
                return Poly_Con_Rel::strictly_intersects();
              if (sp_sign < 0)
                res = Poly_Con_Rel::is_disjoint();
            }
            break;
          } // switch (g.type())
        }
      }
    }
    break;
  } // switch (c.type())
  return res;
}

} // namespace detail

Poly_Con_Rel
Poly_Impl::relation_with(const Con& c) const {
  assert(space_dim() >= c.space_dim());

  if (marked_empty())
    return Poly_Con_Rel::saturates()
      && Poly_Con_Rel::is_included()
      && Poly_Con_Rel::is_disjoint();

  if (space_dim() == 0) {
    if (c.is_inconsistent())
      return (c.is_strict_inequality() && c.inhomo_term() == 0)
        ? Poly_Con_Rel::saturates() && Poly_Con_Rel::is_disjoint()
        : Poly_Con_Rel::is_disjoint();
    else
      return (c.is_equality() || c.inhomo_term() == 0)
        ? Poly_Con_Rel::saturates() && Poly_Con_Rel::is_included()
        : Poly_Con_Rel::is_included();
  }

  ensure_valid_gens();
  if (marked_empty())
    return Poly_Con_Rel::saturates()
      && Poly_Con_Rel::is_included()
      && Poly_Con_Rel::is_disjoint();

  return detail::relation_with(gens(), c);
}

namespace sp {

/* Extend namespace sp by adding methods for Poly_Impl::Sys */

inline bool
satisfied_by_all(const Topol topol, const Con& c,
                 const Poly_Impl::Sys<Gens>& gs) {
  const bool dont_check_ns
    = (topol == Topol::CLOSED)
    || !c.is_strict_inequality()
    || gs.ns_rows.empty();

  if (not satisfied_by_all(c, gs.sing_rows))
    return false;
  if (dont_check_ns)
    return satisfied_by_all(c, gs.sk_rows);

  // Here we have an NNC and a strict inequality;
  // we need to check gs.ns_rows too.
  Bits sat_row;
  if (not compute_sat_row(c, gs.sk_rows, sat_row))
    return false;
  for (const auto& g_ns : gs.ns_rows) {
    if (not g_ns.intersects(sat_row))
      return false;
  }
  // All generators satisfy c.
  return true;
}

inline bool
satisfies_all(const Topol topol, const Gen& g,
              const Poly_Impl::Sys<Cons>& cs) {
  const bool dont_check_ns
    = (topol == Topol::CLOSED) || !g.is_point() || cs.ns_rows.empty();

  if (!satisfies_all(g, cs.sing_rows))
    return false;
  if (dont_check_ns)
    return satisfies_all(g, cs.sk_rows);

  // Here we have an NNC and a point; we need to check cs.ns_rows too.
  Bits sat_row;
  if (!compute_sat_row(g, cs.sk_rows, sat_row))
    return false;
  for (const auto& c_ns : cs.ns_rows) {
    if (!c_ns.intersects(sat_row))
      return false;
  }
  // All constraints are satisfied.
  return true;
}

inline bool
satisfies_all(const Topol topol, const Gen& g,
              const Poly_Impl::Sys<Cons>& cs1,
              const Poly_Impl::Sys<Cons>& cs2) {
  return satisfies_all(topol, g, cs1) && satisfies_all(topol, g, cs2);
}

inline bool
satisfy_all(const Topol topol,
            const Poly_Impl::Sys<Gens>& gs,
            const Poly_Impl::Sys<Cons>& cs) {
  const bool dont_check_ns = (topol == Topol::CLOSED)
    || (gs.ns_rows.empty() && cs.ns_rows.empty());

  // Check gs singular rows (i.e., lines);
  // no need to check them against cs.ns_rows.
  if (!satisfied_by_all(cs.sing_rows, gs.sing_rows))
    return false;
  if (!satisfied_by_all(cs.sk_rows, gs.sing_rows))
    return false;

  // Check gs skel rows.
  if (!satisfied_by_all(cs.sing_rows, gs.sk_rows))
    return false;
  if (dont_check_ns)
    return satisfied_by_all(cs.sk_rows, gs.sk_rows);

  Sat sat(gs.sk_rows.size(), cs.sk_rows.size());
  if (!satisfy_all(gs.sk_rows, cs.sk_rows, sat))
    return false;

  // Check (strict ineqs in) cs.ns_rows ...
  for (const auto& c_ns : cs.ns_rows) {
    // ... against points in gs.sk_rows,
    for (auto i : bwd_index_range(gs.sk_rows)) {
      if (gs.sk_rows[i].is_point() && !c_ns.intersects(sat[i]))
        return false;
    }
    // ... and against points in gs.ns_rows,
    for (const auto& g_ns : gs.ns_rows) {
      if (!c_ns.intersects(detail::sat_all(g_ns, sat)))
        return false;
    }
  }
  // Check strict ineqs in cs.sk_rows against (points in) gs.ns_rows.
  for (auto j : bwd_index_range(cs.sk_rows)) {
    if (cs.sk_rows[j].is_strict_inequality()) {
      Index_Set j_set;
      for (auto i : bwd_index_range(gs.sk_rows)) {
        if (!sat[i][j])
          j_set.set(i);
      }
      if (any_of(gs.ns_rows, [&j_set](const Index_Set& ns) {
                              return subset_eq(ns, j_set);
                             }))
        return false;
    }
  }
  // All checks passed.
  return true;
}

inline bool
satisfy_all(const Topol topol,
            const Poly_Impl::Sys<Gens>& gs1,
            const Poly_Impl::Sys<Gens>& gs2,
            const Poly_Impl::Sys<Cons>& cs1,
            const Poly_Impl::Sys<Cons>& cs2) {
  return satisfy_all(topol, gs1, cs1)
    && satisfy_all(topol, gs1, cs2)
    && satisfy_all(topol, gs2, cs1)
    && satisfy_all(topol, gs2, cs2);
}

} // namespace sp

namespace detail {

inline bool
is_included_in(const Poly_Impl& x, const Poly_Impl& y) {
  assert(x.topol == y.topol);
  assert(x.dim == y.dim);
  assert(!x.marked_empty() && !y.marked_empty());
  assert(x.dim > 0);
  // We need `x' generators.
  x.ensure_valid_gens();
  if (x.marked_empty())
    return true;
  // We need `y' constraints.
  y.ensure_valid_cons();
  // All `x' gens have to satisfy all `y' cons.
  return sp::satisfy_all(x.topol, x.gs, x.gs_pending, y.cs, y.cs_pending);
}

TV_Bool
quick_equals(const Poly_Impl& x, const Poly_Impl& y) {
  assert(x.topol == y.topol);
  assert(x.dim == y.dim);
  assert(!x.marked_empty() && !y.marked_empty());
  assert(x.dim > 0);

  if (!x.marked_min() || !y.marked_min())
    return TV_Bool::DONT_KNOW;

  if (x.num_equals() != y.num_equals()
      || x.num_lines() != y.num_lines()
      || num_rows(x.cs.sk_rows) != num_rows(y.cs.sk_rows)
      || num_rows(x.gs.sk_rows) != num_rows(y.gs.sk_rows)
      || num_rows(x.cs.ns_rows) != num_rows(y.cs.ns_rows)
      || num_rows(x.gs.ns_rows) != num_rows(y.gs.ns_rows))
    return TV_Bool::FALSE;
  /* FIXME.
     In order to enable the following code,
     we need to implement the sorting of the systems
     keeping ns_rows up-to-date.
  if (x.num_equals() == 0) {
    x.cs.sort();
    y.cs.sort();
    return (x.cs == y.cs) ? TV_Bool::TRUE : TV_Bool::FALSE;
  }
  if (x.num_lines() == 0) {
    x.gs.sort();
    y.gs.sort();
    return (x.gs == y.gs) ? TV_Bool::TRUE : TV_Bool::FALSE;
  }
  */
  return TV_Bool::DONT_KNOW;
}

} // namespace detail

Poly_Gen_Rel
Poly_Impl::relation_with(const Gen& g) const {
  assert(space_dim() >= g.space_dim());
  if (is_empty())
    return Poly_Gen_Rel::nothing();
  if (space_dim() == 0)
    return Poly_Gen_Rel::subsumes();
  ensure_valid_cons();
  assert(!marked_empty());
  return sp::satisfies_all(topol, g, cs, cs_pending)
    ? Poly_Gen_Rel::subsumes()
    : Poly_Gen_Rel::nothing();
}

void
Poly_Impl::add_con(Con c) {
  assert(c.space_dim() <= space_dim());
  if (marked_empty())
    return;
  if (space_dim() == 0
      || (c.is_strict_inequality() && is_necessarily_closed())) {
    assert(c.is_tautological() || c.is_inconsistent());
    if (c.is_inconsistent())
      set_empty();
    return;
  }
  ensure_valid_cons();
  assert(marked_min() || has_cs_pending());
  auto& dest = c.is_equality() ? cs_pending.sing_rows : cs_pending.sk_rows;
  dest.push_back(std::move(c));
  if (marked_min())
    set_status(Status::PENDING);
  assert(check_inv());
}

void
Poly_Impl::add_gen(Gen g) {
  assert(g.space_dim() <= space_dim());
  assert(!(is_necessarily_closed() && g.is_closure_point()));
  assert(!(space_dim() == 0 && g.is_line_or_ray()));

  if (space_dim() == 0) {
    if (marked_empty()) {
      assert(g.is_point());
      set_universe();
    }
    assert(check_inv());
    return;
  }
  ensure_valid_gens();

  if (marked_empty()) {
    assert(g.is_point());
    Sys<Gens> tmp;
    tmp.sk_rows.push_back(std::move(g));
    reinit_with_gens(dim, topol, std::move(tmp));
    return;
  }
  assert(marked_min() || has_gs_pending());
  auto& dest = g.is_line() ?  gs_pending.sing_rows : gs_pending.sk_rows;
  dest.push_back(std::move(g));
  if (marked_min())
    set_status(Status::PENDING);
  assert(check_inv());
}

namespace detail {

// Returns pairs <i1, i2>, meaning that rows[i1] makes rows[i2] a duplicate.
template <typename Rows>
inline std::vector<std::pair<dim_type, dim_type>>
find_duplicates(const Rows& rows) {
  const auto sz = num_rows(rows);
  assert(sz > 0);
  Dims sorted(sz);
  std::iota(sorted.begin(), sorted.end(), 0);
  std::sort(sorted.begin(), sorted.end(),
            [&rows](dim_type i1, dim_type i2) {
              return compare(rows[i1], rows[i2]) < 0;
            });
  std::vector<std::pair<dim_type, dim_type>> dupl;
  for (auto i1 = 0; i1 != sz - 1; ++i1) {
    for (auto i2 = i1 + 1; i2 != sz; ++i2) {
      if (compare(rows[i1], rows[i2]) == 0)
        dupl.emplace_back(i1, i2);
      else {
        i1 = i2 - 1;
        break;
      }
    }
  }
  return dupl;
}

template <typename Sys>
inline void
maybe_filter_duplicates(Sys& sys, dim_type threshold) {
  if (threshold == not_a_dim() || num_rows(sys.sk_rows) <= threshold)
    return;
  auto sk_dupl = find_duplicates(sys.sk_rows);
  if (sk_dupl.empty())
    return;
  Index_Set tbr;
  for (const auto& d : sk_dupl)
    tbr.set(d.second);
  erase_using_sorted_indices(sys.sk_rows, tbr);
  // Rewrite the non-skeleton elements.
  for (auto& ns : sys.ns_rows) {
    if (!ns.intersects(tbr))
      continue;
    for (const auto& d : sk_dupl)
      if (ns.test(d.second)) {
        ns.set(d.first);
        ns.reset(d.second);
      }
  }
  // Now sort and unique non-skel elements.
  std::sort(sys.ns_rows.begin(), sys.ns_rows.end());
  sys.ns_rows.erase(std::unique(sys.ns_rows.begin(), sys.ns_rows.end()),
                    sys.ns_rows.end());
}

} // namespace detail

void
Poly_Impl::minimize() const {
  if (marked_empty() || marked_min() || (space_dim() == 0))
    return;

  // The polyhedron is only `logically' const.
  auto& x = const_cast<Poly_Impl&>(*this);

  if (has_valid_cons()) {
    detail::maybe_filter_duplicates(x.cs_pending, minimize_filter_threshold);
    if (sat_c.num_rows() == 0)
      x.sat_c = x.sat_g.transpose();
    if (detail::add_and_minimize<true>(x.dim, x.cs, x.cs_pending,
                                       x.gs, x.sat_c, x.sat_g))
      x.set_empty();
    else
      x.set_status(Status::MINIMIZED);
    // check_inv() would trigger an infinite recursion.
    assert(check_inv_without_minimize());
    return;
  }

  assert(has_valid_gens());
  detail::maybe_filter_duplicates(x.gs_pending, minimize_filter_threshold);
  if (sat_g.num_rows() == 0)
    x.sat_g = x.sat_c.transpose();
  if (detail::add_and_minimize<false>(x.dim, x.gs, x.gs_pending,
                                      x.cs, x.sat_g, x.sat_c))
    PPLITE_UNREACH;
  else
    x.set_status(Status::MINIMIZED);
  // check_inv() would trigger an infinite recursion.
  assert(check_inv_without_minimize());
}

namespace detail {

// If the transformation is not invertible we may transform
// valid lines and rays into the origin of the space.
// Fixing up the system is up to the caller.
void
affine_image_unsafe(Gens& gs, Var var, const Linear_Expr& expr,
                    const Integer& inhomo, const Integer& den) {
  assert(den > 0);
  Integer num;
  const auto expr_nz = expr.non_zeroes();
  for (auto& g : gs) {
    auto& g_expr = g.impl().expr;
    auto& g_inhomo = g.impl().inhomo;
    sp::assign(num, expr_nz, expr, inhomo, g_expr, g_inhomo);
    if (den != 1) {
      g_expr *= den;
      g_inhomo *= den;
    }
    g_expr.set(var, num);
    g.strong_normalize();
  }
}

void
affine_image(Poly_Impl::Sys<Gens>& gs, Var var,
             Linear_Expr expr, Integer inhomo, Integer den) {
  if (den < 0) {
    neg_assign(expr);
    neg_assign(inhomo);
    neg_assign(den);
  }
  affine_image_unsafe(gs.sing_rows, var, expr, inhomo, den);
  affine_image_unsafe(gs.sk_rows, var, expr, inhomo, den);
  const bool invertible = (expr.get(var) != 0);
  if (!invertible)
    remove_invalid_lines_and_rays(gs);
}

void
affine_preimage(Cons& cs, Var var, const Linear_Expr& expr,
                const Integer& inhomo, const Integer& den) {
  assert(den > 0);
  const auto& var_e = expr.get(var);
  for (auto& c : cs) {
    if (c.coeff(var) == 0)
      continue;
    auto var_c = c.coeff(var);
    auto& c_expr = c.impl().expr;
    auto& c_inhomo = c.impl().inhomo;
    if (den != 1) {
      c_expr *= den;
      c_inhomo *= den;
    }
    // FIXME: CHECKME.
    add_mul_assign(c_expr, var_c, expr);
    add_mul_assign(c_inhomo, var_c, inhomo);
    c_expr.set(var, var_c * var_e);
    c.strong_normalize();
  }
}

void
affine_preimage(Poly_Impl::Sys<Cons>& cs, Var var,
                Linear_Expr expr, Integer inhomo, Integer den) {
  if (den < 0) {
    neg_assign(expr);
    neg_assign(inhomo);
    neg_assign(den);
  }
  affine_preimage(cs.sing_rows, var, expr, inhomo, den);
  affine_preimage(cs.sk_rows, var, expr, inhomo, den);
}

} // namespace detail

void
Poly_Impl::affine_image(const Var var, const Linear_Expr& expr,
                        const Integer& inhomo, const Integer& den) {
  assert(space_dim() > 0);
  assert(var.space_dim() <= space_dim());
  assert(expr.space_dim() <= space_dim());
  assert(den != 0);

  if (marked_empty())
    return;
  const bool invertible = (expr.get(var) != 0);
  if (invertible) {
    detail::affine_image(gs, var, expr, inhomo, den);
    if (!gs_pending.empty())
      detail::affine_image(gs_pending, var, expr, inhomo, den);
    // Build the inverse transformation.
    const auto& inv_den = expr.get(var);
    auto inv_expr = -expr;
    inv_expr.set(var, den);
    auto inv_inhomo = -inhomo;
    detail::affine_preimage(cs, var, inv_expr, inv_inhomo, inv_den);
    if (!cs_pending.empty())
      detail::affine_preimage(cs_pending, var, inv_expr, inv_inhomo, inv_den);
    assert(check_inv());
    return;
  }
  // The transformation is not invertible.
  assert(!invertible);
  unconstrain(var);
  add_con(den * var - expr == inhomo);
}

void
Poly_Impl::affine_preimage(const Var var, const Linear_Expr& expr,
                           const Integer& inhomo, const Integer& den) {
  assert(space_dim() > 0);
  assert(var.space_dim() <= space_dim());
  assert(expr.space_dim() <= space_dim());
  assert(den != 0);
  if (marked_empty())
    return;
  const bool invertible = (expr.get(var) != 0);
  if (invertible) {
    // Invertible: compute inverse and delegate to affine_image.
    const auto& inv_den = expr.get(var);
    auto inv_expr = -expr;
    inv_expr.set(var, den);
    auto inv_inhomo = -inhomo;
    affine_image(var, inv_expr, inv_inhomo, inv_den);
    return;
  }
  // The transformation is not invertible.
  assert(!invertible);
  add_con(den * var - expr == inhomo);
  unconstrain(var);
}

namespace detail {

bool
check_par_affine(dim_type space_dim,
                 const Vars& vars,
                 const Linear_Exprs& exprs,
                 const Integers& inhomos, const Integers& dens) {
  const auto n_vars = vars.size();
  if (n_vars == 0)
    return false;
  if (n_vars != exprs.size()
      || n_vars != inhomos.size()
      || n_vars != dens.size())
    return false;

  std::vector<bool> is_target(space_dim, false);
  for (auto j : bwd_index_range(vars)) {
    if (vars[j].space_dim() > space_dim)
      return false;
    if (exprs[j].space_dim() > space_dim)
      return false;
    if (dens[j] == 0)
      return false;
    const auto v_id = vars[j].id();
    if (is_target[v_id])
      return false;
    is_target[v_id] = true;
  }
  return true;
}

} // namespace detail

void
Poly_Impl::parallel_affine_image(const Vars& vars,
                                 const Linear_Exprs& exprs,
                                 const Integers& inhomos,
                                 const Integers& dens) {
  assert(detail::check_par_affine(dim, vars, exprs, inhomos, dens));
  detail::par_affine_image_aux(*this, vars, exprs, inhomos, dens);
  assert(check_inv());
}

namespace detail {

// Remove strict inequalities preserving the ray-face-cutter.
void
safe_remove_strict_ineqs(Cons& cs_sk, NS_Rows& cs_ns,
                         const Gens& gs_sk, const Sat& sat_c) {
  // Remove skeleton strict ineqs, except for positivity, if present.
  strict_ineqs_become_nonstrict_ineqs(cs_sk);
  // Remove non-skeleton strict ineqs
  cs_ns.clear();

  auto pos_i = std::find_if(cs_sk.begin(), cs_sk.end(),
                            std::mem_fn(&Con::is_tautological));
  const bool has_pos = (pos_i != cs_sk.end());
  if (has_pos) {
    // Positivity stays strict.
    assert(has_strict_ineq(cs_sk));
    return;
  }

  // Here positivity is not present in cs_sk.
  // We must reintroduce the rfc in cs_ns.
  Index_Set rays;
  for (auto g : bwd_index_range(gs_sk))
    if (gs_sk[g].is_ray())
      rays.set(g);
  Bits rays_sat = detail::sat_all(rays, sat_c);
  rays_sat.complement_until(num_rows(cs_sk));
  Index_Set efc(std::move(rays_sat));
  cs_ns.push_back(std::move(efc));
}

} // namespace detail

void
Poly_Impl::topological_closure_assign() {
  if (is_necessarily_closed() || marked_empty() || (space_dim() == 0))
    return;
  ensure_valid_gens();
  if (marked_empty())
    return;
  gs.ns_rows.clear();
  detail::closure_points_become_points(gs.sk_rows);
  gs_pending.ns_rows.clear();
  detail::closure_points_become_points(gs_pending.sk_rows);
  assert(!has_cs_pending());
  if (marked_min() || has_gs_pending()) {
    detail::safe_remove_strict_ineqs(cs.sk_rows, cs.ns_rows,
                                     gs.sk_rows, sat_c);
  }
  assert(check_inv());
}

void
Poly_Impl::intersection_assign(const Poly_Impl& y) {
  auto& x = *this;
  assert(x.topol == y.topol);
  assert(x.dim == y.dim);
  if (x.marked_empty())
    return;
  if (y.marked_empty()) {
    x.set_empty();
    return;
  }
  // Handle 0-dim universe case.
  if (x.dim == 0)
    return;

  x.ensure_valid_cons();
  assert(!x.marked_empty());
  y.minimize();
  assert(y.check_inv());
  if (y.marked_empty()) {
    x.set_empty();
    return;
  }
  assert(x.marked_min() || x.has_cs_pending());
  assert(y.marked_min());
  detail::concat_sys(x.cs_pending, y.cs);

  x.set_status(Status::PENDING);
  assert(x.check_inv());
}

void
Poly_Impl::poly_hull_assign(const Poly_Impl& y) {
  auto& x = *this;
  assert(x.topol == y.topol);
  assert(x.dim == y.dim);
  if (x.marked_empty()) {
    x = y;
    return;
  }
  if (x.dim == 0 || y.marked_empty())
    return;

  x.ensure_valid_gens();
  if (x.marked_empty()) {
    x = y;
    return;
  }
  y.minimize();
  assert(y.check_inv());
  if (y.marked_empty())
    return;
  assert(x.marked_min() || x.has_gs_pending());
  assert(y.marked_min());
  detail::concat_sys(x.gs_pending, y.gs);

  x.set_status(Status::PENDING);
  assert(x.check_inv());
}

namespace detail {

void
concatenate_minimized(Poly_Impl& x, const Poly_Impl& y) {
  assert(x.space_dim() > 0 && y.space_dim() > 0);
  assert(!x.is_empty() && !y.is_empty());
  assert(x.is_minimized() && y.is_minimized());

  const bool not_closed
    = !x.is_topologically_closed() || !y.is_topologically_closed();

  // Gather info on x and y empty face cutters
  auto x_info = get_pos_efc_info(x);
  const auto x_pos_index = x_info.first;
  const auto x_efc_index = x_info.second;
  auto y_info = get_pos_efc_info(y);
  const auto y_pos_index = y_info.first;
  const auto y_efc_index = y_info.second;

  const bool drop_x_pos = (x_pos_index != not_a_dim())
    && (y_pos_index == not_a_dim());
  const bool drop_x_efc = (x_efc_index != not_a_dim())
    && (y_pos_index == not_a_dim())
    && (y_efc_index == not_a_dim());
  const bool drop_y_efc = (x_pos_index == not_a_dim());

  // Get rid of x pos/efc if it will become redundant.
  if (drop_x_pos) {
    Index_Set tbr {x_pos_index};
    detail::remove_rows(tbr, x.cs.sk_rows, x.cs.ns_rows, x.sat_g);
    x.sat_c = x.sat_g.transpose();
  } else if (drop_x_efc)
    x.cs.ns_rows.erase(x.cs.ns_rows.begin() + x_efc_index);

  // Save size info (after purging pos, but before concatenating).
  const auto x_cs_num_skel = num_rows(x.cs.sk_rows);
  const auto y_cs_num_skel = num_rows(y.cs.sk_rows);

  // Concatenate constraint systems.
  {
    // Concatenate singular cons.
    x.cs.sing_rows.reserve(x.cs.sing_rows.size() + y.cs.sing_rows.size());
    for (const auto& y_sg : y.cs.sing_rows) {
      x.cs.sing_rows.emplace_back(y_sg);
      x.cs.sing_rows.back().shift_space_dims(0, x.dim);
    }
    // Concatenate skel cons (except y pos).
    x.cs.sk_rows.reserve(x_cs_num_skel + y_cs_num_skel);
    for (auto i : range(y_cs_num_skel)) {
      const auto& y_sk = y.cs.sk_rows[i];
      if (i == y_pos_index) {
        assert(y_sk.is_strict_inequality() && y_sk.is_tautological());
        continue;
      }
      x.cs.sk_rows.emplace_back(y_sk);
      x.cs.sk_rows.back().shift_space_dims(0, x.dim);
    }
    // Concatenate non-skel cons.
    x.cs.ns_rows.reserve(x.cs.ns_rows.size() + y.cs.ns_rows.size());
    for (auto i : index_range(y.cs.ns_rows)) {
      const auto& y_ns = y.cs.ns_rows[i];
      if (i == y_efc_index && drop_y_efc)
        continue;
      x.cs.ns_rows.emplace_back(y_ns);
      x.cs.ns_rows.back() >>= x_cs_num_skel;
    }
    // Maybe adjust x efc to match concatenation.
    if (x_efc_index != not_a_dim() && y_efc_index != not_a_dim()) {
      auto new_efc = y.cs.ns_rows[y_efc_index];
      new_efc >>= x_cs_num_skel;
      auto& x_efc = x.cs.ns_rows[x_efc_index];
      x_efc |= new_efc;
    }
  } // Concatenate constraint systems


  { // Concatenate generator systems.
    // Helper: compute sets of rays and (closure) points.
    auto rays_and_cps = [](const Gens& gs,
                           Index_Set& rays, Index_Set& points, Index_Set& cps) {
      for (auto i : bwd_index_range(gs)) {
        assert(!gs[i].is_line());
        if (gs[i].is_ray())
          rays.set(i);
        else {
          cps.set(i);
          if (gs[i].is_point())
            points.set(i);
        }
      }
    };

    // Helper: concat g1 with g2 (the latter to be shifted by sdim1)
    auto gen_concat = [](Gen& g1, const Gen& g2, dim_type sdim1) {
      assert(g1.is_point_or_closure_point() && g2.is_point_or_closure_point());
      if (g1.is_point() && g2.is_closure_point())
        g1.set_type(Gen::CLOSURE_POINT);
      auto& e1 = g1.linear_expr();
      const auto& e2 = g2.linear_expr();
      Integer lcm;
      lcm_assign(lcm, g1.divisor(), g2.divisor());
      const auto mult1 = lcm / g1.divisor();
      const auto mult2 = lcm / g2.divisor();
      e1 *= mult1;
      g1.impl().inhomo *= mult1;
      for (auto i : bwd_dim_range(e2)) {
        if (e2.get(i).is_zero())
          continue;
        add_mul_assign(e1, mult2 * e2.get(i), Var(sdim1 + i));
      }
      g1.strong_normalize();
    };

    // Helper: concat r1 with r2; the bits of r2 are copied after
    // the size1 bits of r1, discarding (if present) the bit corresponding
    // to the skel positivity constraint in r2 (having index `skip2').
    auto sat_concat = [](const Bits& r1, dim_type size1,
                         const Bits& r2, dim_type skip2) {
      assert(size1 >= 0 && (skip2 == not_a_dim() || skip2 >= 0));
      if (r2.empty())
        return r1;
      Bits res = r2;
      if (skip2 == not_a_dim() || r2.last() < skip2) {
        // No need to erase bit skip2
        res >>= size1;
        if (!r1.empty())
          res.copy_from(0, r1, 0, r1.last() + 1);
        return res;
      }
      // Here we need to erase skip2
      // FIXME: find better way
      res.reset(skip2);
      if (res.empty())
        return r1;
      res >>= (size1 - 1);
      if (size1 > 0)
        res.reset(size1 - 1);
      if (skip2 > 0)
        res.copy_from(size1, r2, 0, skip2);
      if (!r1.empty())
        res.copy_from(0, r1, 0, r1.last() + 1);
      return res;
    };

    // Concatenate singular gens (lines).
    x.gs.sing_rows.reserve(x.gs.sing_rows.size() + y.gs.sing_rows.size());
    for (const auto& y_sg : y.gs.sing_rows) {
      x.gs.sing_rows.emplace_back(y_sg);
      x.gs.sing_rows.back().shift_space_dims(0, x.dim);
    }

    // Concatenate skel gens (producing new_skel and new_sat).
    Gens new_sk;
    Sat new_sat;
    const auto& x_gs_sk = x.gs.sk_rows;
    const auto& y_gs_sk = y.gs.sk_rows;

    // These will be used to remap ns gens (if not closed)
    using XY_Map
      = std::unordered_map<std::pair<dim_type, dim_type>,
                           dim_type, Pair_Hasher<dim_type, dim_type>>;
    XY_Map xy_map;
    std::unordered_map<dim_type, dim_type> xr_map, yr_map;

    // Compute the sets of rays and (closure) points.
    Index_Set x_rays, x_points, x_cps;
    rays_and_cps(x_gs_sk, x_rays, x_points, x_cps);
    Index_Set y_rays, y_points, y_cps;
    rays_and_cps(y_gs_sk, y_rays, y_points, y_cps);

    // Avoid reallocations
    const auto new_rows = x_rays.size() + y_rays.size()
      + (x_cps.size() * y_cps.size());
    const auto new_cols = x.cs.sk_rows.size();
    new_sk.reserve(new_rows);
    new_sat.impl().rows.reserve(new_rows);
    new_sat.add_cols(new_cols);

    // Concatenate the (closure) points in x_gs_sk with those in y_gs_sk.
    for (auto i : x_cps) {
      const auto& x_i = x_gs_sk[i];
      assert(!x_i.is_ray());
      for (auto j : y_cps) {
        if (not_closed)
          xy_map[std::make_pair(i,j)] = new_sk.size();
        const auto& y_j = y_gs_sk[j];
        assert(!y_j.is_ray());
        new_sk.emplace_back(x_i);
        gen_concat(new_sk.back(), y_j, x.dim);
        new_sat.add_row(sat_concat(x.sat_c[i], x.sat_c.num_cols(),
                                   y.sat_c[j], y_pos_index));
      }
    }
    // Add x rays (can move them)
    for (auto i : x_rays) {
      if (not_closed)
        xr_map[i] = new_sk.size();
      auto& x_i = x.gs.sk_rows[i];
      assert(x_i.is_ray());
      new_sk.emplace_back(std::move(x_i));
      new_sat.add_row(std::move(x.sat_c[i]));
    }
    // Add (and shift) y rays
    for (auto j : y_rays) {
      if (not_closed)
        yr_map[j] = new_sk.size();
      const auto& y_j = y_gs_sk[j];
      assert(y_j.is_ray());
      new_sk.emplace_back(y_j);
      new_sk.back().shift_space_dims(0, x.dim);
      new_sat.add_row(sat_concat(Bits(), x.sat_c.num_cols(),
                                 y.sat_c[j], y_pos_index));
    }

    // (Maybe) Concatenate non-skel gens (producing new_ns).
    NS_Rows new_ns;
    if (not_closed) {
      // Helper function: concatenation with index remapping.
      auto ns_concat = [&xy_map, &xr_map, &yr_map, &x_rays, &y_rays]
        (const Index_Set& x_ns, const Index_Set& y_ns) {
        Index_Set res;
        for (auto i : x_ns) {
          if (x_rays.test(i)) {
            assert(xr_map.find(i) != xr_map.end());
            res.set(xr_map[i]);
            continue;
          }
          for (auto j : y_ns) {
            if (y_rays.test(j)) {
              assert(yr_map.find(j) != xr_map.end());
              res.set(yr_map[j]);
              continue;
            }
            assert(xy_map.find(std::make_pair(i,j)) != xy_map.end());
            res.set(xy_map[std::make_pair(i,j)]);
          }
        }
        return res;
      };

      // Avoid reallocations.
      new_ns.reserve((x_points.size() + x.gs.ns_rows.size())
                     * (y_points.size() + y.gs.ns_rows.size()));
      // Combine x ns points with y ns points
      for (const auto& x_ns : x.gs.ns_rows)
        for (const auto& y_ns : y.gs.ns_rows)
          new_ns.emplace_back(ns_concat(x_ns, y_ns));
      // Combine x skel points with y ns points
      for (auto i : x_points) {
        auto x_ns = Index_Set(i);
        for (const auto& y_ns : y.gs.ns_rows)
          new_ns.emplace_back(ns_concat(x_ns, y_ns));
      }
      // Combine y skel points with x ns points
      for (auto j : y_points) {
        auto y_ns = Index_Set(j);
        for (const auto& x_ns : x.gs.ns_rows)
          new_ns.emplace_back(ns_concat(x_ns, y_ns));
      }
    }

    // Move new sk/ns gens and sat_c matrix in their place.
    x.gs.sk_rows = std::move(new_sk);
    x.gs.ns_rows = std::move(new_ns);
    x.sat_c = std::move(new_sat);
  } // Concatenate generators

  // Update space dimension and sat_g.
  x.dim += y.dim;
  x.sat_g = x.sat_c.transpose();
  assert(x.check_inv());
}

} // namespace detail

void
Poly_Impl::concatenate_assign(const Poly_Impl& y) {
  auto& x = *this;
  assert(x.topol == y.topol);
  if (x.is_empty()) {
    x.dim += y.dim;
    return;
  }
  if (y.is_empty()) {
    x.dim += y.dim;
    x.set_empty();
    return;
  }
  if (y.space_dim() == 0)
    return;
  if (x.space_dim() == 0) {
    x = y;
    return;
  }
  x.minimize();
  y.minimize();
  detail::concatenate_minimized(x, y);
}

void
Poly_Impl::time_elapse_assign(const Poly_Impl& y) {
  Poly_Impl& x = *this;
  assert(x.topol == y.topol && x.dim == y.dim);
  // Deal with the zero-dim special case.
  if (x.dim == 0) {
    if (y.marked_empty())
      x.set_empty();
    return;
  }
  // Deal with the empty special case.
  if (x.marked_empty() || y.marked_empty()
      || x.is_empty() || y.is_empty()) {
    x.set_empty();
    return;
  }

  assert(x.has_valid_gens() && y.has_valid_gens());
  auto& x_lines = x.gs_pending.sing_rows;
  std::copy(y.gs.sing_rows.begin(), y.gs.sing_rows.end(),
            std::back_inserter(x_lines));
  std::copy(y.gs_pending.sing_rows.begin(), y.gs_pending.sing_rows.end(),
            std::back_inserter(x_lines));
  auto& x_sk = x.gs_pending.sk_rows;
  detail::add_as_rays(y.gs.sk_rows, x_sk);
  detail::add_as_rays(y.gs_pending.sk_rows, x_sk);
  // Note: y non-skel points can be ignored.
  if (x.marked_min() && !x.gs_pending.empty())
    x.set_status(Status::PENDING);
  assert(x.check_inv());
}

void
Poly_Impl::widening_assign(const Poly_Impl& y, const Cons* upto_ptr,
                           Widen_Impl w_impl, Widen_Spec w_spec) {
  auto& x = *this;
  assert(x.topol == y.topol);
  assert(x.dim == y.dim);

  // We work on minimized representation.
  x.minimize();
  y.minimize();
  if (detail::widening_preamble(x, y, w_spec))
    return;

  // Before modifying x, select valid upto constraints.
  Index_Set valid = detail::valid_upto_cons(x, upto_ptr);
  if (w_spec == Widen_Spec::SAFE) {
    // In the safe case, we have to check y too.
    valid &= detail::valid_upto_cons(y, upto_ptr);
  }
  // Apply specific widening algorithm.
  switch (w_impl) {
  case Widen_Impl::H79:
    h79_widen::h79_widen(x, y, w_spec);
    break;
  case Widen_Impl::BOXED_H79:
    h79_widen::boxed_h79_widen(x, y, w_spec);
    break;
  case Widen_Impl::BHRZ03:
    bhrz03_widen::bhrz03_widen(x, y, w_spec);
    break;
  }
  // Now add valid upto constraints.
  detail::add_valid_upto_cons(x, valid, upto_ptr);
}

void
Poly_Impl::add_space_dims(dim_type m, bool project) {
  if (m == 0)
    return;
  if (marked_empty()) {
    dim += m;
    return;
  }
  if (project) {
    /* Projecting on zero. */
    ensure_valid_cons();
    auto& eqs = cs.sing_rows;
    eqs.reserve(eqs.size() + m);
    for (auto i : range(m))
      eqs.push_back(Var(dim + i) == 0);
    dim += m;
    assert(check_inv());
    return;
  }
  /* Embedding. */
  ensure_valid_gens();
  if (marked_empty()) {
    dim += m;
    return;
  }
  auto& lines = gs.sing_rows;
  lines.reserve(dim + m);
  for (auto i : range(m))
    lines.push_back(line(Var(dim + i)));
  dim += m;
  assert(check_inv());
}

void
Poly_Impl::remove_higher_space_dims(dim_type new_dim) {
  assert(new_dim <= space_dim());
  if (new_dim == space_dim())
    return;
  Dims tbr(space_dim() - new_dim);
  std::iota(tbr.begin(), tbr.end(), new_dim);
  remove_space_dims(tbr.begin(), tbr.end());
}

namespace detail {

Linear_Expr
map_space_dims(Linear_Expr&& src, const Dims& pfunc, dim_type new_dim) {
  Linear_Expr res;
  res.set_space_dim(new_dim);
  const auto i_end = std::min(src.space_dim(), num_rows(pfunc));
  for (auto i : range(i_end)) {
    if (pfunc[i] != not_a_dim())
      res[pfunc[i]] = std::move(src[i]);
  }
  return res;
}

void
map_space_dims_unsafe(Gens& gens, const Dims& pfunc, dim_type new_dim) {
  // Note: this function may create invalid rays/lines.
  // We do not eliminate them here (it is done by the caller),
  // because we need to keep the non-skel component consistent.
  for (auto& g : gens) {
    auto& e = g.linear_expr();
    e = map_space_dims(std::move(e), pfunc, new_dim);
    g.strong_normalize();
  }
}

void
map_space_dims(Poly_Impl::Sys<Gens>& sys,
               const Dims& pfunc, dim_type new_dim) {
  map_space_dims_unsafe(sys.sing_rows, pfunc, new_dim);
  map_space_dims_unsafe(sys.sk_rows, pfunc, new_dim);
  remove_invalid_lines_and_rays(sys);
}

} // namespace detail

void
Poly_Impl::map_space_dims(const Dims& pfunc) {
  assert(check_inv());
  assert(space_dim() == num_rows(pfunc));
  if (dim == 0)
    return;
  const auto max_dim = *std::max_element(pfunc.begin(), pfunc.end());
  if (max_dim == not_a_dim()) {
    const bool nonempty = !is_empty();
    dim = 0;
    if (nonempty)
      set_universe();
    assert(check_inv());
    return;
  }

  const auto new_dim = 1 + max_dim;
  if (is_empty()) {
    dim = new_dim;
    assert(check_inv());
    return;
  }

  if (new_dim == dim) {
    // pfunc is a permutation
    detail::permute_space_dims(*this, pfunc, dim);
    assert(check_inv());
    return;
  }

  // `pfunc' is not a permutation.
  // Implement it as remove_space_dims + permutation.
  Dims tbr;
  Dims perm;
  for (auto i : range(dim)) {
    if (pfunc[i] == not_a_dim())
      tbr.push_back(i);
    else
      perm.push_back(pfunc[i]);
  }
  remove_space_dims(tbr.begin(), tbr.end());
  detail::permute_space_dims(*this, perm, dim);
  assert(check_inv());
}

namespace detail {

template <typename Rows>
bool
check_inv(const Poly_Impl::Sys<Rows>& sys, dim_type dim) {
  using Row = typename Rows::value_type;
  if (any_of(sys.sing_rows, [dim](const Row& sg) {
                             return !is_singular(sg)
                              || (sg.space_dim() > dim)
                              || !(sg.check_inv());
                            }))
    return false;
  if (any_of(sys.sk_rows, [dim](const Row& sk) {
                           return is_singular(sk)
                            || (sk.space_dim() > dim)
                            || !(sk.check_inv());
                          }))
    return false;
  const auto sk_size = num_rows(sys.sk_rows);
  if (any_of(sys.ns_rows, [sk_size](const Index_Set& ns) {
                           return (ns.size() < 2) || (ns.last() >= sk_size);
                          }))
    return false;
  // All checks passed.
  return true;
}

/**
 * There must always be a cutter.
 * - If there are rays, there must be the ray-face-cutter,
 *   i.e., a constraint rfc such that: sat_g(rfc) = R.
 *   The strict positivity constraint, when present, is a valid rfc.
 *   When it is not present, it is implied by some constraints in cs_sk;
 *   these constraints are the rfc.
 * - If there is no ray, the efc is the maximal one.
 */
bool
check_efc(const NS_Rows& cs_ns,
          const Cons& cs_sk,
          const Gens& gs_sk,
          const Sat& sat_g, const Sat& sat_c) {
  if (!detail::check_cutters(cs_ns, cs_sk, gs_sk, sat_g))
    return false;

  auto efc_i = std::find_if(cs_ns.cbegin(), cs_ns.cend(),
                            [&sat_g, &gs_sk](const Index_Set& ns) {
                              return is_empty_face_cutter(ns, sat_g, gs_sk);
                            });
  // We already know that there is a valid cutter.
  if (efc_i == cs_ns.cend())
    return true;

  // Check that the efc is the correct rfc.
  // Compute the expected rfc.
  Index_Set rays;
  for (auto g : bwd_index_range(gs_sk))
    if (gs_sk[g].is_ray())
      rays.set(g);
  Bits rfc_sat = sat_all(rays, sat_c);
  rfc_sat.complement_until(num_rows(cs_sk));
  // Note: if `rays' is empty, this computes the maximal set.
  Index_Set real_efc(std::move(rfc_sat));

  return (real_efc == *efc_i);
}

// Find the empty-face-cutter
// either in maximal or rfc representation.
// If it is made redundant returns an empty Index_Set.
Index_Set
find_efc(const NS_Rows& cs_ns,
         const Cons& cs_sk, const Gens& gs_sk,
         const Sat& sat_g, const Sat& sat_c) {
  Index_Set empty_efc;
  // Compute the set of rays.
  Index_Set rays;
  for (dim_type i : bwd_index_range(gs_sk))
    if (gs_sk[i].is_ray())
      rays.set(i);
  auto maybe_efc_sat = detail::sat_all(rays, sat_c);
  maybe_efc_sat.complement_until(num_rows(cs_sk));
  const Index_Set maybe_efc(std::move(maybe_efc_sat));

  const auto efc_size = maybe_efc.size();
  const auto cs_sk_size = num_rows(cs_sk);
  if (efc_size == 1) {
    // It is the positivity constraint (singleton).
    assert(cs_sk[*maybe_efc.begin()].is_strict_inequality());
    assert(cs_sk[*maybe_efc.begin()].is_tautological());

    // Check that the found set cuts the empty face.
    assert(is_empty_face_cutter(maybe_efc, sat_g, gs_sk));
    return maybe_efc;
  }

  if (efc_size == cs_sk_size) {
    // It is maximal.
    // It can be made redundant by any other strict constraint.
    if (cs_ns.size() == 0 || cs_ns[0].size() != cs_sk_size) {
      assert(cs_ns.size() > 0 ||
             any_of(cs_sk, std::mem_fn(&Con::is_strict_inequality)));
      return empty_efc;
    }

    // Check that the found set cuts the empty face.
    assert(maybe_efc == cs_ns[0]);
    assert(rays.empty());
    assert(is_empty_face_cutter(maybe_efc, sat_g, gs_sk));
    return maybe_efc;
  }

  // It is the rays face cutter.
  assert(!rays.empty());
  auto has_strict = [&cs_sk](dim_type i)
    { return cs_sk[i].is_strict_inequality(); };
  auto has_ns = [&maybe_efc](const Index_Set& ns)
    { return subset_ne(ns, maybe_efc); };
  // The rfc can be made redundant by some strict skel constraint.
  if (any_of(maybe_efc, has_strict))
    return empty_efc;
  // The rfc can be made redundant by some non-skel constraint.
  if (any_of(cs_ns, has_ns))
    return empty_efc;

  // The rfc is found. Check that is present in ns_rows.
  assert(any_of(cs_ns, [&maybe_efc](const Index_Set& ns)
                       { return ns == maybe_efc; }));
  // Check that the found set cuts the empty face.
  assert(is_empty_face_cutter(maybe_efc, sat_g, gs_sk));
  // Silence warning (sat_g is used only when in debugging mode).
  (void) sat_g;
  return maybe_efc;
}

void
maybe_dump(const Poly_Impl& ph, const char* reason) {
#ifdef NDEBUG
  (void) ph;
  (void) reason;
#else // In debugging mode, be noisy.
  std::cerr << reason << std::endl;
  ph.ascii_dump(std::cerr);
  std::cerr << reason << std::endl;
#endif
}

} // namespace detail

bool
Poly_Impl::check_inv() const {
  const char* reason = aux_check_inv_without_minimize();
  if (reason != nullptr) {
    detail::maybe_dump(*this, reason);
    return false;
  }
  reason = aux_check_inv_check_minimize();
  if (reason != nullptr) {
    detail::maybe_dump(*this, reason);
    return false;
  }
  return true;
}

bool
Poly_Impl::check_inv_without_minimize() const {
  const char* reason = aux_check_inv_without_minimize();
  if (reason != nullptr) {
    detail::maybe_dump(*this, reason);
    return false;
  }
  return true;
}

const char*
Poly_Impl::aux_check_inv_without_minimize() const {
  /* Check 0-dim polyhedra. */
  if (dim == 0 && not marked_empty() && not marked_min())
    return "Zero-dim poly neither empty nor minimized";

  /* Check systems cardinalities. */
  if (not detail::check_inv(cs, dim))
    return "non-pending constraints are broken";
  if (not detail::check_inv(gs, dim))
    return "non-pending generators are broken";
  if (not detail::check_inv(cs_pending, dim))
    return "pending constraints are broken";
  if (not detail::check_inv(gs_pending, dim))
    return "pending generators are broken";

  // Checks status consistency.
  switch (status) {

  case Status::EMPTY:
    if (not ((cs.sing_rows.size() == 1)
             && cs.sing_rows[0].space_dim() == 0
             && cs.sk_rows.empty() && cs.ns_rows.empty()
             && cs_pending.empty()
             && gs.empty() && gs_pending.empty()
             && sat_c.empty() && sat_g.empty()))
      return "poly marked empty is broken";
    break;

  case Status::PENDING:
    if (cs_pending.empty() && gs_pending.empty())
      return "poly marked pending has no pending cons/gens";
    if (not cs_pending.empty() && not gs_pending.empty())
      return "poly marked pending has both pending cons and gens";
    if (cs.empty())
      return "poly marked pending has empty non-pending cons";
    if (gs.empty())
      return "poly marked pending has empty non-pending gens";
    // Intentionally falling through:
    // we also check minimization invariants for the non-pending part.

    [[ fallthrough ]];
  case Status::MINIMIZED:
    if (not (num_rows(cs.sing_rows) <= dim
             && num_rows(gs.sing_rows) <= dim
             && not sat_c.empty() && not sat_g.empty()
             && (sat_c.num_rows() == sat_g.num_cols())
             && (sat_g.num_rows() == sat_c.num_cols())
             && (sat_c.num_rows() == num_rows(gs.sk_rows))
             && (sat_g.num_rows() == num_rows(cs.sk_rows))))
      return "minimized poly: broken cardinalities";
    // Note: recheck status (may be falling through from PENDING).
    if (status == Status::MINIMIZED
        && not (cs_pending.empty() && gs_pending.empty()))
      return "minimized poly has pending cons/gens";
    // Check sat matrix consistency.
    if (not (sat_g == sat_c.transpose()))
      return "poly: sat_c is not the transpose of sat_g";
    if (not detail::check_sat(sat_g, cs.sk_rows, gs.sk_rows))
      return "poly: sat matrix does not match scalar prod";
    // Check for efc
    if (not detail::check_efc(cs.ns_rows, cs.sk_rows,
                              gs.sk_rows, sat_g, sat_c))
      return "poly: invalid efc";
    break;

  default:
    return "poly: unknown status value";
  } // switch (status)

  // All checks passed
  return nullptr;
}

const char*
Poly_Impl::aux_check_inv_check_minimize() const {
  /* This (computationally heavy) check assumes that
     aux_check_inv_without_minimize() has succeded.
     Since it calls the following Poly_Impl methods:
       1) empty constructor
       2) universe constructor
       3) copy constructor
       4) reinit_with_gens
       5) minimize
       6) equals
     Methods 2-6 cannot call this check, neither directly nor indirectly
     by calling check_inv(), since this would trigger an infinite recursion.
  */
  if (status != Status::MINIMIZED && status != Status::PENDING)
    return nullptr;

  auto get_counters = [](const Poly_Impl& ph) {
    auto num_sk_strict
      = std::count_if(ph.cs.sk_rows.begin(), ph.cs.sk_rows.end(),
                      std::mem_fn(&Con::is_strict_inequality));
    auto num_sk_points
      = std::count_if(ph.gs.sk_rows.begin(), ph.gs.sk_rows.end(),
                      std::mem_fn(&Gen::is_point));
    auto num_ns_strict = num_rows(ph.cs.ns_rows);
    auto num_ns_points = num_rows(ph.gs.ns_rows);
    return std::make_tuple(num_sk_strict, num_sk_points,
                           num_ns_strict, num_ns_points);
  };

  // Rebuild non-pending part from constraints.
  {
    // Copy the non-pending part of this into x.
    Poly_Impl x = *this;
    x.cs_pending.clear();
    x.gs_pending.clear();
    x.status = Status::MINIMIZED;
    // Create y using the constraints of x as pending.
    Poly_Impl y(x.dim, Spec_Elem::UNIVERSE, x.topol);
    y.cs_pending = x.cs;
    y.set_status(Status::PENDING);
    y.minimize();

    const auto x_counters = get_counters(x);
    const auto y_counters = get_counters(y);
    if (x_counters != y_counters) {
      return "poly said minimized, but it is not"
        ": conversion cons->gens computes a DD pair having"
        " a cardinality mismatch for cons/gens kinds";
    }
    if (not y.equals(x)) {
      return "poly said minimized, but it is not"
        ": conversion cons->gens computes different polyhedron";
    }
  }
  // Rebuild non-pending part from generators.
  {
    // Copy the non-pending part of this into x.
    Poly_Impl x = *this;
    x.cs_pending.clear();
    x.gs_pending.clear();
    x.status = Status::MINIMIZED;
    // Create y using the generators of x as pending.
    Poly_Impl y(x.dim, Spec_Elem::UNIVERSE, x.topol);
    auto x_gs_copy = x.gs;
    y.reinit_with_gens(x.dim, x.topol, std::move(x_gs_copy));
    y.minimize();

    const auto x_counters = get_counters(x);
    const auto y_counters = get_counters(y);
    if (x_counters != y_counters) {
      return "poly said minimized, but it is not"
        ": conversion gens->cons computes a DD pair having"
        " a cardinality mismatch for cons/gens kinds";
    }
    if (not y.equals(x)) {
      return "poly said minimized, but it is not"
        ": conversion gens->cons computes different polyhedron";
    }
  }
  // All checks passed
  return nullptr;
}

void
Poly_Impl::ascii_dump(std::ostream& s) const {
  using namespace IO_Operators;
  s << "topol " << topol << "\n";
  s << "dim " << dim << "\n";
  s << "status " << status << "\n";
  s << "=> cs sys\n";
  cs.ascii_dump(s);
  s << "=> gs sys\n";
  gs.ascii_dump(s);
  s << "sat_c\n";
  sat_c.ascii_dump(s);
  s << "sat_g\n";
  sat_g.ascii_dump(s);
  s << "=> cs_pending\n";
  cs_pending.ascii_dump(s);
  s << "=> gs_pending\n";
  gs_pending.ascii_dump(s);
}

namespace detail {

template <typename Row>
bool
ascii_load(std::vector<Row>& rows, dim_type nrows, std::istream& s) {
  rows.resize(nrows);
  for (auto i : index_range(rows)) {
    if (!rows[i].ascii_load(s))
      return false;
  }
  return true;
}

template <typename Rows>
bool
ascii_load(Poly_Impl::Sys<Rows>& sys, std::istream& s) {
  sys.clear();
  dim_type num = 0;
  if (not (ascii_load_string(s, "sing_rows")
           && (s >> num) && (num >= 0)
           && ascii_load(sys.sing_rows, num, s)))
    return false;
  if (not (ascii_load_string(s, "sk_rows")
           && (s >> num) && (num >= 0)
           && ascii_load(sys.sk_rows, num, s)))
    return false;
  if (not (ascii_load_string(s, "ns_rows")
           && (s >> num) && (num >= 0)
           && ascii_load(sys.ns_rows, num, s)))
    return false;
  return true;
}

} // namespace detail

bool
Poly_Impl::ascii_load(std::istream& s) {
  std::string str;
  if (not (ascii_load_string(s, "topol") && (s >> str)))
    return false;
  if (str == "CLOSED")
    topol = Topol::CLOSED;
  else if (str == "NNC")
    topol = Topol::NNC;
  else
    return false;

  if (not (ascii_load_string(s, "dim") && (s >> dim) && (dim >= 0)))
    return false;

  if (not (ascii_load_string(s, "status") && (s >> str)))
    return false;
  if (str == "EMPTY")
    status = Status::EMPTY;
  else if (str == "MINIMIZED")
    status = Status::MINIMIZED;
  else if (str == "PENDING")
    status = Status::PENDING;
  else
    return false;

  if (not (ascii_load_string(s, "=>")
           && ascii_load_string(s, "cs")
           && ascii_load_string(s, "sys")
           && detail::ascii_load(cs, s)))
    return false;
  if (not (ascii_load_string(s, "=>")
           && ascii_load_string(s, "gs")
           && ascii_load_string(s, "sys")
           && detail::ascii_load(gs, s)))
    return false;

  if (not (ascii_load_string(s, "sat_c") && sat_c.ascii_load(s)))
    return false;
  if (not (ascii_load_string(s, "sat_g") && sat_g.ascii_load(s)))
    return false;

  if (not (ascii_load_string(s, "=>")
           && ascii_load_string(s, "cs_pending")
           && detail::ascii_load(cs_pending, s)))
    return false;
  if (not (ascii_load_string(s, "=>")
           && ascii_load_string(s, "gs_pending")
           && detail::ascii_load(gs_pending, s)))
    return false;

  assert(check_inv());
  return true;
}

void
Poly_Impl::print(std::ostream& s) const {
  if (is_empty()) {
    s << "false";
    return;
  }
  minimize();
  bool comma = false;
  for (const auto& c : cons()) {
    if (comma) s << ", ";
    using namespace IO_Operators;
    s << c;
    comma = true;
  }
}

namespace IO_Operators {

std::ostream&
operator<<(std::ostream& s, Topol topol) {
  s << (topol == Topol::CLOSED ? "C" : "NNC");
  return s;
}

std::ostream&
operator<<(std::ostream& s, Poly_Impl::Status status) {
  using Status = Poly_Impl::Status;
  switch (status) {
  case Status::EMPTY:
    s << "EMPTY";
    break;
  case Status::MINIMIZED:
    s << "MINIMIZED";
    break;
  case Status::PENDING:
    s << "PENDING";
    break;
  default:
    PPLITE_UNREACH;
    break;
  }
  return s;
}

std::ostream&
operator<<(std::ostream& s, const Poly_Impl& ph) {
  ph.print(s);
  return s;
}

std::ostream&
operator<<(std::ostream& s, const Poly& ph) {
  s << ph.impl();
  return s;
}

} // namespace IO_Operators

namespace detail {

// Just for debugging purposes (pass by value is meant).
bool quick_equals_test(Poly_Impl x, Poly_Impl y) {
  x.minimize();
  y.minimize();
  // quick_equals assumes polyhedra are not 0-dim and are not empty.
  if (x.dim == 0 || x.is_empty() || y.is_empty())
    return true;
  auto res = quick_equals(x, y);
  if (res == TV_Bool::DONT_KNOW)
    return true;
  assert(res != TV_Bool::TRUE);
  // If quick_equals says false, it has to be correct.
  return (!is_included_in(y, x) || !is_included_in(x, y));
}

} // namespace detail

bool
Poly_Impl::equals(const Poly_Impl& y) const {
  const auto& x = *this;
  if (x.topol != y.topol || x.dim != y.dim)
    return false;

  if (x.marked_empty())
    return y.is_empty();
  else if (y.marked_empty())
    return x.is_empty();
  else if (x.space_dim() == 0)
    return true;

  assert(detail::quick_equals_test(x, y));
  switch (detail::quick_equals(x, y)) {
  case TV_Bool::TRUE:
    return true;
  case TV_Bool::FALSE:
    return false;
  case TV_Bool::DONT_KNOW:
  default:
    using detail::is_included_in;
    if (!is_included_in(x, y))
      return false;
    if (x.marked_empty())
      return y.is_empty();
    else
      return is_included_in(y, x);
  }
}

bool
Poly_Impl::contains(const Poly_Impl& y) const {
  const auto& x = *this;
  assert(x.topol == y.topol);
  assert(x.dim == y.dim);
  if (y.marked_empty())
    return true;
  if (x.marked_empty())
    return y.is_empty();
  if (y.space_dim() == 0)
    return true;
  if (detail::quick_equals(x, y) == TV_Bool::TRUE)
    return true;
  return detail::is_included_in(y, x);
}

bool
Poly_Impl::is_disjoint_from(const Poly_Impl& y) const {
  const auto& x = *this;
  if (x.marked_empty() || y.marked_empty())
    return true;
  if (x.is_empty() || y.is_empty())
    return true;

  // Speculative optimizations.
  auto nx = x.num_min_cons();
  if (nx < 2)
    return detail::simple_poly_is_disjoint_from_poly(nx, x, y);
  auto ny = y.num_min_cons();
  if (ny < 2)
    return detail::simple_poly_is_disjoint_from_poly(ny, y, x);

  // General case: compute intersection (trying to exploit incrementality).
  if (nx >= ny) {
    Poly_Impl z = x;
    z.intersection_assign(y);
    return z.is_empty();
  } else {
    Poly_Impl z = y;
    z.intersection_assign(x);
    return z.is_empty();
  }
}

namespace detail {

// THIS IS THE OLD VERSION, KEPT JUST FOR DEBUGGING PURPOSES
BBox get_bounding_box(const Poly_Impl& ph) {
  if (ph.marked_empty())
    return BBox(ph.space_dim(), Spec_Elem::EMPTY);
  if (ph.space_dim() == 0)
    return BBox(ph.space_dim(), Spec_Elem::UNIVERSE);

  ph.minimize();
  if (ph.is_empty())
    return BBox(ph.space_dim(), Spec_Elem::EMPTY);

  auto res = BBox(ph.space_dim());
  res.set_origin();
  // Scan lines.
  for (const auto& g : ph.gs.sing_rows)
    detail::add_as_line(res, g);
  // Scan rays.
  for (auto r_i : bwd_index_range(ph.gs.sk_rows)) {
    const auto& r = ph.gs.sk_rows[r_i];
    if (r.is_ray())
      detail::add_as_ray(res, r);
  }
  // Now scan skeleton (closure) points.
  bool point_seen = false;
  for (auto g_i : bwd_index_range(ph.gs.sk_rows)) {
    const auto& g = ph.gs.sk_rows[g_i];
    if (!g.is_point() && !g.is_closure_point())
      continue;
    const auto& div = g.divisor();
    if (point_seen) {
      // This is not the first (closure) point.
      detail::add_as_point(res, g);
    } else {
      // This is the first (closure) point seen.
      point_seen = true;
      // Note: do not use detail::init_with_point(),
      // as here we want to preserve infinite bounds.
      for (auto i : bwd_dim_range(ph)) {
        if (not res.constrains(Var(i)))
          continue;
        Rational value(g.coeff(Var(i)), div);
        auto& itv = res[i];
        if (itv.has_lb())
          itv.set_lb(value);
        if (itv.has_ub())
          itv.set_ub(std::move(value));
      }
    }
  }
  res.maybe_update_volume_info();
  assert(res.check_inv());
  return res;
}

} // namespace detail

BBox
Poly_Impl::get_bounding_box() const {
  if (marked_empty())
    return BBox(space_dim(), Spec_Elem::EMPTY);
  if (space_dim() == 0)
    return BBox(space_dim(), Spec_Elem::UNIVERSE);

  minimize();
  if (is_empty())
    return BBox(space_dim(), Spec_Elem::EMPTY);

  auto res = BBox(space_dim());
  res.set_origin();
  Index_Set done_lb;
  Index_Set done_ub;
  // Search for interval constraints.
  for (const auto& c : cs.sing_rows) {
    if (!is_proper_interval_con(c))
      continue;
    const auto& ex = c.linear_expr();
    auto i = ex.last_nonzero();
    res[i].set_singleton(Rational(-c.inhomo_term(), ex.get(i)));
    done_lb.set(i);
    done_ub.set(i);
  }
  for (const auto& c : cs.sk_rows) {
    if (!is_proper_interval_con(c))
      continue;
    const auto& ex = c.linear_expr();
    auto i = ex.last_nonzero();
    const auto& coeff = ex.get(i);
    if (sgn(coeff) > 0) {
      res[i].lb = Rational(-c.inhomo_term(), coeff);
      done_lb.set(i);
    } else {
      assert(sgn(coeff) < 0);
      res[i].ub = Rational(-c.inhomo_term(), coeff);
      done_ub.set(i);
    }
  }
  // Scan lines.
  for (const auto& g : gs.sing_rows) {
    for (auto i : bwd_dim_range(g))
      if (g.coeff(Var(i)) != 0) {
        res[i].set_universe();
        done_lb.set(i);
        done_ub.set(i);
      }
  }
  // Scan rays.
  for (auto r_i : bwd_index_range(gs.sk_rows)) {
    const auto& r = gs.sk_rows[r_i];
    if (!r.is_ray())
      continue;
    for (auto i : dim_range(r)) {
        auto s = sgn(r.coeff(Var(i)));
        if (s < 0) {
          res[i].unset_lb();
          done_lb.set(i);
        }
        else if (s > 0) {
          res[i].unset_ub();
          done_ub.set(i);
        }
      }
  }
  // Now scan skeleton (closure) points.
  for (auto i : bwd_dim_range(*this)) {
    if (done_lb.test(i) && done_ub.test(i))
      continue;
    bool point_seen = false;
    for (auto g_i : bwd_index_range(gs.sk_rows)) {
      const auto& g = gs.sk_rows[g_i];
      if (!g.is_point() && !g.is_closure_point())
        continue;
      const auto& div = g.divisor();
      Rational value(g.coeff(Var(i)), div);
      auto& itv = res[i];
      if (point_seen) {
        // This is not the first (closure) point.
        // Update values
        if (!done_lb.test(i) && value < itv.lb)
          itv.set_lb(value);
        if (!done_ub.test(i) && value > itv.ub)
          itv.set_ub(std::move(value));
      } else {
        // This is the first (closure) point seen.
        point_seen = true;
        if (!done_lb.test(i))
          itv.set_lb(value);
        if (!done_ub.test(i))
          itv.set_ub(std::move(value));
      }
    } // for each g
  } // for each dim
  res.maybe_update_volume_info();
  assert(res.check_inv());
  assert(res == detail::get_bounding_box(*this));
  return res;
}

bool
Poly_Impl::boxed_contains(const Poly_Impl& y) const {
  const auto& x = *this;
  assert(x.topol == y.topol);
  assert(x.dim == y.dim);
#ifndef NDEBUG
  {
    // Precondition: the topologically closed
    // bounding box of x contains that of y.
    auto x_box = x.get_bounding_box();
    auto y_box = y.get_bounding_box();
    assert(x_box.contains(y_box));
  }
#endif
  x.minimize();
  y.minimize();
  if (x.dim == 0 || x.marked_empty() || y.marked_empty())
    return true;

  // The check can be limited to the constraints of x
  // that are NOT closed interval constraints.

  // Singular constraints (i.e., equalities)
  for (const auto& c : x.cs.sing_rows) {
    if (is_proper_interval_con(c)) {
      assert(sp::satisfied_by_all(x.topol, c, y.gs));
      continue;
    }
    if (!sp::satisfied_by_all(x.topol, c, y.gs))
      return false;
  }
  // Skel constraints.
  for (const auto& c : x.cs.sk_rows) {
    if (c.is_nonstrict_inequality() && is_proper_interval_con(c)) {
      assert(sp::satisfied_by_all(x.topol, c, y.gs));
      continue;
    }
    if (!sp::satisfied_by_all(x.topol, c, y.gs))
      return false;
  }
  // Non-skel constraints.
  if (x.is_necessarily_closed())
    return true;
  // Avoid useless materialization of strict positivity
  // (only when it is maximal, which is often the case).
  bool has_maximal_efc
    = x.cs.ns_rows.size() == 1
    && x.cs.ns_rows[0].size() == num_rows(x.cs.sk_rows);
  if (has_maximal_efc)
    return true;
  // Materialize non-skel constraints.
  for (const auto& ns : x.cs.ns_rows) {
    Con c = detail::materialize(ns, x.cs.sk_rows);
    if (!sp::satisfied_by_all(x.topol, c, y.gs))
      return false;
  }
  return true;
}

bool
Poly_Impl::is_bounded_expr(bool from_below, const Linear_Expr& expr) const {
  // A zero-dim or empty polyhedron makes everything bounded.
  if (space_dim() == 0 || is_empty())
    return true;

  // Speculative optimization: usually only a few dims are non-zero.
  const auto expr_nz = expr.non_zeroes();

  // Predicate checking if a line/ray causes expr to be unbounded.
  auto makes_unbounded = [&expr_nz, &expr, from_below](const Gen& g) {
    auto sp_sign = sp::sign(expr_nz, expr, g.linear_expr());
    if (sp_sign == 0)
      return false;
    if (g.is_line())
      return true;
    assert(g.is_ray());
    return not static_cast<bool>(from_below ^ (sp_sign < 0));
  };

  // Check lines.
  if (any_of(gs.sing_rows, makes_unbounded))
    return false;
  if (any_of(gs_pending.sing_rows, makes_unbounded))
    return false;
  // Check rays.
  for (const auto& g : gs.sk_rows) {
    if (g.is_ray() && makes_unbounded(g))
      return false;
  }
  for (const auto& g : gs_pending.sk_rows) {
    if (g.is_ray() && makes_unbounded(g))
      return false;
  }
  // expr is bounded.
  return true;
}

bool
Poly_Impl::constrains(Var var) const {
  assert(var.space_dim() <= space_dim());
  if (marked_empty())
    return true;

  auto c_pred = [var](const Con& c) {
    return c.linear_expr().get(var) != 0;
  };
  auto cs_pred = [c_pred](const Cons& cs) {
    return any_of(cs, c_pred);
  };
  if (marked_min())
    return cs_pred(cs.sing_rows) || cs_pred(cs.sk_rows);

  if (has_valid_gens()) {
    // Try a quick, incomplete check for the universe polyhedron.
    if (num_rows(gs.sing_rows) == space_dim())
      return false;
    // Predicate for lines.
    auto lines_pred = [var](const Gen& g) {
      assert(g.is_line());
      return g.linear_expr().all_zeroes_except(var);
    };
    // Scan lines.
    if (any_of(gs.sing_rows, lines_pred))
      return false;
    // Scan pending lines.
    if (any_of(gs_pending.sing_rows, lines_pred))
      return false;

    // Scan skel (also pending) for opposite rays.
    bool has_pos_ray = false;
    bool has_neg_ray = false;
    // Predicate for rays (note: it has status).
    auto rays_pred = [var, &has_pos_ray, &has_neg_ray](const Gen& g) {
      if (g.is_ray() && g.linear_expr().all_zeroes_except(var)) {
        const int sign = sgn(g.coeff(var));
        if (sign > 0) {
          has_pos_ray = true;
          if (has_neg_ray)
            return true;
        } else {
          assert(sign < 0);
          has_neg_ray = true;
          if (has_pos_ray)
            return true;
        }
      }
      return false;
    };
    // Scan skel for opposite rays.
    if (any_of(gs.sk_rows, rays_pred))
      return false;
    // Scan pending skel for opposite rays.
    if (any_of(gs_pending.sk_rows, rays_pred))
      return false;
  }

  minimize();
  if (marked_empty())
    return true;
  return cs_pred(cs.sing_rows) || cs_pred(cs.sk_rows);
}

Index_Set
Poly_Impl::get_unconstrained() const {
  Index_Set res;
  if (is_empty())
    return res;
  minimize();
  for (const auto& line : gs.sing_rows) {
    const auto& le = line.linear_expr();
    if (is_proper_interval_expr(le))
      res.set(le.first_nonzero());
  }
  return res;
}

bool
Poly_Impl::min(const Affine_Expr& ae, Rational& value,
               bool* included_ptr, Gen* g_ptr) const {
  assert(ae.space_dim() <= space_dim());
  if (space_dim() == 0) {
    if (marked_empty())
      return false;
    value = Rational(ae.inhomo);
    if (included_ptr)
      *included_ptr = true;
    if (g_ptr)
      *g_ptr = point();
    return true;
  }
  if (is_empty())
    return false;

  assert(has_valid_gens());
  auto ae_nz = ae.expr.non_zeroes();

  // Check lines: returns true if g makes expr (lower) unbounded.
  auto check_line = [&ae_nz, &ae](const Gen& g) {
    auto sp_sign = sp::sign(ae_nz, ae.expr, g.linear_expr());
    return (sp_sign != 0);
  };

  if (any_of(gs.sing_rows, check_line))
    return false;
  if (any_of(gs_pending.sing_rows, check_line))
    return false;

  // Check skel and non-skel generators.
  bool need_init = true;
  const Gen* cand_ptr = nullptr;

  // Helper to see if we need to also check for equal-value generators
  // (this implies the need to check for non-skel points).
  auto need_extra_checks = [&]() {
    return (included_ptr && not *included_ptr)
      || (cand_ptr && cand_ptr->is_closure_point());
  };

  // Check skel: returns true if g makes expr (lower) unbounded;
  // updates current minimum info otherwise.
  auto check_skel = [&](const Gen& g) {
    Integer sp;
    sp::add_sp_assign(sp, ae_nz, ae.expr, g.linear_expr());
    if (g.is_ray()) {
      if (sgn(sp) < 0)
        return true;
    } else {
      assert(g.is_point_or_closure_point());
      Rational new_value(sp, g.divisor());
      if (need_init || new_value < value) {
        value = std::move(new_value);
        if (included_ptr)
          *included_ptr = g.is_point();
        if (g_ptr)
          cand_ptr = &g;
        need_init = false;
      } else if (need_extra_checks()
                 && g.is_point()
                 && (new_value == value)) {
        // Need anyway to update included_ptr and/or cand.
        if (included_ptr)
          *included_ptr = true;
        if (g_ptr)
          cand_ptr = &g;
      }
    }
    return false;
  };

  // Skel component.
  if (any_of(gs.sk_rows, check_skel))
    return false;
  if (any_of(gs_pending.sk_rows, check_skel))
    return false;

  if (not need_extra_checks()) {
    value += Rational(ae.inhomo);
    if (g_ptr) {
      assert(cand_ptr);
      *g_ptr = *cand_ptr;
    }
    return true;
  }

  // Non-skel component: materialize.
  Gen mater_cand;
  for (const auto& ns : gs.ns_rows) {
    Gen mater = detail::materialize(ns, gs.sk_rows);
    if (check_skel(mater))
      PPLITE_UNREACH;
    if (cand_ptr == &mater) {
      mater_cand = std::move(mater);
      cand_ptr = &mater_cand;
    }
  }
  for (const auto& ns : gs_pending.ns_rows) {
    Gen mater = detail::materialize(ns, gs_pending.sk_rows);
    if (check_skel(mater))
      PPLITE_UNREACH;
    if (cand_ptr == &mater) {
      mater_cand = std::move(mater);
      cand_ptr = &mater_cand;
    }
  }

  value += Rational(ae.inhomo);
  if (g_ptr) {
    assert(cand_ptr);
    *g_ptr = *cand_ptr;
  }
  return true;
}

Itv
Poly_Impl::get_bounds(Var var) const {
  assert(var.space_dim() <= space_dim());
  if (is_empty())
    return Itv::empty();
  assert(has_valid_gens());

  // Check lines.
  for (const auto& g : gs.sing_rows) {
    if (not g.coeff(var).is_zero())
      return Itv::universe();
  }
  // Check pending lines.
  for (const auto& g : gs_pending.sing_rows) {
    if (not g.coeff(var).is_zero())
      return Itv::universe();
  }

  bool inf_lb = false;
  bool inf_ub = false;
  // Returns true if g makes the result universe
  auto check_ray = [var, &inf_lb, &inf_ub](const Gen& g) {
    if (not g.is_ray())
      return false;
    auto sg = sgn(g.coeff(var));
    if (sg < 0) {
      inf_lb = true;
      return inf_ub;
    } else if (sg > 0) {
      inf_ub = true;
      return inf_lb;
    }
    return false;
  };

  // Check rays.
  for (const auto& g : gs.sk_rows) {
    check_ray(g);
    if (inf_lb && inf_ub)
      return Itv::universe();
  }
  // Check pending rays.
  for (const auto& g : gs_pending.sk_rows) {
    check_ray(g);
    if (inf_lb && inf_ub)
      return Itv::universe();
  }

  bool first = true;
  Itv res;
  auto check_point = [var, inf_lb, inf_ub, &first, &res](const Gen& g) {
    if (g.is_ray())
      return;
    Rational value = Rational(g.coeff(var), g.divisor());
    if (first) {
      if (inf_lb) {
        assert(!inf_ub);
        res.set_ub(std::move(value));
      } else if (inf_ub) {
        assert(!inf_lb);
        res.set_lb(std::move(value));
      } else
        res.set_singleton(std::move(value));
      first = false;
    } else {
      if (!inf_lb && value < res.lb)
        res.set_lb(std::move(value));
      else if (!inf_ub && value > res.ub)
        res.set_ub(std::move(value));
    }
  };

  // Check skel (closure) points.
  for (const auto& g : gs.sk_rows)
    check_point(g);
  // Check pending skel (closure) points.
  for (const auto& g : gs_pending.sk_rows)
    check_point(g);
  // No need to check non-skel points.
  return res;
}

Itv
Poly_Impl::get_bounds(const Affine_Expr& ae) const {
  assert(ae.space_dim() <= space_dim());
  if (is_empty())
    return Itv::empty();
  assert(has_valid_gens());

  const auto ae_nz = ae.expr.non_zeroes();

  // Returns true if g makes the result universe
  auto check_line = [&ae_nz, &ae](const Gen& g) {
    assert(g.is_line());
    auto sp_sign = sp::sign(ae_nz, ae.expr, g.linear_expr());
    return (sp_sign != 0);
  };

  // Check lines.
  for (const auto& g : gs.sing_rows) {
    if (check_line(g))
      return Itv::universe();
  }
  // Check pending lines.
  for (const auto& g : gs_pending.sing_rows) {
    if (check_line(g))
      return Itv::universe();
  }

  bool inf_lb = false;
  bool inf_ub = false;
  // Returns true if g makes the result universe
  auto check_ray = [&ae_nz, &ae, &inf_lb, &inf_ub](const Gen& g) {
    assert(not (inf_lb && inf_ub));
    if (not g.is_ray())
      return false;
    auto sp_sign = sp::sign(ae_nz, ae.expr, g.linear_expr());
    if (sp_sign == 0)
      return false;
    else if (sp_sign == 1) {
      inf_ub = true;
      if (inf_lb)
        return true;
    } else {
      inf_lb = true;
      if (inf_ub)
        return true;
    }
    return false;
  };

  // Check rays.
  for (const auto& g : gs.sk_rows) {
    if (check_ray(g))
      return Itv::universe();
  }
  // Check pending rays.
  for (const auto& g : gs_pending.sk_rows) {
    if (check_ray(g))
      return Itv::universe();
  }

  bool first = true;
  Itv res;
  // Note: this cannot make res universe.
  auto check_point
    = [inf_lb, inf_ub, &ae_nz, &ae, &first, &res](const Gen& g) {
    if (g.is_ray())
      return;
    Integer num;
    sp::add_sp_assign(num, ae_nz, ae.expr, g.linear_expr());
    auto value = Rational(std::move(num), g.divisor());
    if (first) {
      if (inf_lb) {
        assert(!inf_ub);
        res.set_ub(std::move(value));
      } else if (inf_ub) {
        assert(!inf_lb);
        res.set_lb(std::move(value));
      } else
        res.set_singleton(std::move(value));
      first = false;
    } else {
      if (!inf_lb && value < res.lb)
        res.set_lb(std::move(value));
      else if (!inf_ub && value > res.ub)
        res.set_ub(std::move(value));
    }
  };

  // Check skel (closure) points.
  for (const auto& g : gs.sk_rows)
    check_point(g);
  // Check pending skel (closure) points.
  for (const auto& g : gs_pending.sk_rows)
    check_point(g);
  // No need to check non-skel points.
  // Finally, take into account inhomo term.
  if (not ae.inhomo.is_zero()) {
    auto r = Rational(ae.inhomo);
    if (res.has_lb())
      res.lb += r;
    if (res.has_ub())
      res.ub += r;
  }
  return res;
}

Itv
Poly_Impl::get_bounds(const Itv_Expr& ie) const {
  if (is_empty())
    return Itv::empty();

  // Computes scalar product on itvs and expr.
  auto scalar_prod = [](const Itv_Expr& ie, const Linear_Expr& expr) {
    const auto& ie_vars = ie.first;
    const auto& ie_itvs = ie.second;
    assert(num_rows(ie_vars) == num_rows(ie_itvs));

    auto res = Itv::zero();
    for (auto i : index_range(ie_vars)) {
      auto ie_dim = ie_vars[i].id();
      const auto& ie_itv = ie_itvs[i];
      const auto& c = expr.get(ie_dim);
      if (c.is_zero() || ie_itv.is_zero())
        continue;
      if (c.is_one())
        res.add_assign(ie_itv);
      else {
        Itv prod = ie_itv;
        prod.mul_assign(Rational(c));
        res.add_assign(prod);
      }
    }
    return res;
  };

  // Returns true if g makes the result universe
  auto check_line = [&scalar_prod, &ie](const Gen& g) {
    assert(g.is_line());
    auto sp = scalar_prod(ie, g.linear_expr());
    return not sp.is_zero();
  };

  // Check lines.
  for (const auto& g : gs.sing_rows) {
    if (check_line(g))
      return Itv::universe();
  }
  // Check pending lines.
  for (const auto& g : gs_pending.sing_rows) {
    if (check_line(g))
      return Itv::universe();
  }

  bool inf_lb = false;
  bool inf_ub = false;
  // Returns true if g makes the result universe
  auto check_ray = [&scalar_prod, &ie, &inf_lb, &inf_ub](const Gen& g) {
    assert(not (inf_lb && inf_ub));
    if (not g.is_ray())
      return false;
    auto sp = scalar_prod(ie, g.linear_expr());
    assert(sp.is_bounded());
    if (sp.is_zero())
      return false;
    if (sgn(sp.ub) == 1) {
      inf_ub = true;
      if (inf_lb)
        return true;
    }
    if (sgn(sp.lb) == -1) {
      inf_lb = true;
      if (inf_ub)
        return true;
    }
    return false;
  };

  // Check rays.
  for (const auto& g : gs.sk_rows) {
    if (check_ray(g))
      return Itv::universe();
  }
  // Check pending rays.
  for (const auto& g : gs_pending.sk_rows) {
    if (check_ray(g))
      return Itv::universe();
  }

  bool first = true;
  Itv res;
  // Note: this cannot make res universe.
  auto check_point
    = [&scalar_prod, inf_lb, inf_ub, &ie, &first, &res](const Gen& g) {
    if (g.is_ray())
      return;
    auto sp = scalar_prod(ie, g.linear_expr());
    if (not g.divisor().is_one())
      sp.mul_assign(Rational(1, g.divisor()));
    if (not first) {
      res.lub_assign(sp);
    } else {
      first = false;
      res = std::move(sp);
      if (inf_lb) {
        assert(!inf_ub);
        res.unset_lb();
      } else if (inf_ub) {
        assert(!inf_lb);
        res.unset_ub();
      }
    }
  };

  // Check skel (closure) points.
  for (const auto& g : gs.sk_rows)
    check_point(g);
  // Check pending skel (closure) points.
  for (const auto& g : gs_pending.sk_rows)
    check_point(g);
  // No need to check non-skel points.
  return res;
}

namespace detail {

Index_Set
fold_gens(Gens& gens, const Index_Set& vars, Var v) {
  Gens new_gens;
  Index_Set invalid_gens;
  auto is_valid = [](const Gen& g) {
    return g.is_point_or_closure_point() || !g.linear_expr().is_zero();
  };
  const auto old_dst = v.id();
  const auto new_dst = old_dst
    - std::count_if(vars.begin(), vars.end(),
                    [old_dst](dim_type id) { return id < old_dst; });
  const auto min_dim = new_dst + 1;

  Index_Set vv = vars;
  vv.set(v.id());
  for (auto k : index_range(gens)) {
    auto& g = gens[k];
    // Compute min and max coeffs for dimensions in vars.
    auto mm = std::minmax_element(vv.begin(), vv.end(),
                                  [&g](dim_type i, dim_type j) {
                                    return g.coeff(Var(i)) < g.coeff(Var(j));
                                  });
    // Note: copies are meant.
    Integer c_min = g.coeff(Var(*mm.first));
    Integer c_max = g.coeff(Var(*mm.second));
    g.linear_expr().remove_space_dims(vars.begin(), vars.end());
    if (g.space_dim() < min_dim)
      g.set_space_dim(min_dim);
    if (c_min != c_max) {
      // Note: copy has to be done *before* modifying g.
      Gen new_g = g;
      new_g.linear_expr()[new_dst] = std::move(c_max);
      if (is_valid(new_g)) {
        new_g.strong_normalize();
        new_gens.push_back(std::move(new_g));
      }
    }
    g.linear_expr()[new_dst] = std::move(c_min);
    if (is_valid(g))
      g.strong_normalize();
    else
      invalid_gens.set(k);
  }
  // Move new_gens into gens (at the end).
  gens.insert(gens.end(), new_gens.begin(), new_gens.end());
  return invalid_gens;
}

} // namespace detail

void
Poly_Impl::fold_space_dims(const Index_Set& vars, Var dst_var) {
  assert(dst_var.space_dim() <= space_dim());
  assert(!vars.test(dst_var.id()));
  if (vars.empty())
    return;
  assert(vars.last() < space_dim());

  // Generators are needed.
  ensure_valid_gens();
  // Update space dimension.
  dim -= vars.size();
  if (marked_empty())
    return;
  if (dim == 0) {
    set_universe();
    return;
  }

  // FIXME: speculative optimization.
  // If status == MINIMIZED, then we could check if the equalities
  // imply that all vars are equal to dst_var.
  // In this case, folding is a nop.

  // Move everything into new gs system.
  Sys<Gens> tmp = std::move(gs);
  detail::concat_sys(tmp, std::move(gs_pending));
  // Materialize non-skel generators.
  // CHECKME: is there any more efficient way?
  for (const auto& ns : tmp.ns_rows)
    tmp.sk_rows.push_back(detail::materialize(ns, tmp.sk_rows));
  tmp.ns_rows.clear();

  // Fold dims in geometric generators.
  Index_Set invalid = detail::fold_gens(tmp.sing_rows, vars, dst_var);
  erase_using_sorted_indices(tmp.sing_rows, invalid);
  invalid = detail::fold_gens(tmp.sk_rows, vars, dst_var);
  erase_using_sorted_indices(tmp.sk_rows, invalid);
  reinit_with_gens(dim, topol, std::move(tmp));
}

void
Poly_Impl::expand_space_dim(Var var, dim_type m) {
  if (m == 0)
    return;
  if (is_empty()) {
    dim += m;
    return;
  }

  // FIXME: optimizable.
  Cons new_cons;
  for (const auto& c : cons()) {
    if (c.coeff(var).is_zero())
      continue;
    for (auto d : range(m)) {
      Con new_c = c;
      new_c.linear_expr().swap_space_dims(var.id(), dim + d);
      new_cons.push_back(std::move(new_c));
    }
  }

  add_space_dims(m);
  add_cons(new_cons);
  assert(check_inv());
}

Poly_Impl
Poly_Impl::split(const Con& c, Topol t) {
  assert(t == Topol::CLOSED || t == topology());
  assert(!c.is_equality());
  assert(t == Topol::NNC || c.is_nonstrict_inequality());
  auto& ph_then = *this;
  ph_then.minimize();
  auto ph_else = Poly_Impl(ph_then.space_dim(),
                           Spec_Elem::EMPTY,
                           ph_then.topology());
  if (not ph_then.is_empty())
    detail::Q_split_poly(ph_then, ph_else, c, t);
  assert(ph_then.check_inv() && ph_else.check_inv());
  assert(ph_then.is_minimized() && ph_else.is_minimized());
  return ph_else;
}

Poly_Impl
Poly_Impl::integral_split(const Con& c) {
  assert(is_necessarily_closed());
  auto& ph_then = *this;
  ph_then.minimize();
  auto ph_else = Poly_Impl(ph_then.space_dim(),
                           Spec_Elem::EMPTY,
                           ph_then.topology());
  if (not ph_then.is_empty())
    detail::Z_split_poly(ph_then, ph_else, c);
  assert(ph_then.check_inv() && ph_else.check_inv());
  return ph_else;
}

// Initialize static (thread local) data member.
PPLITE_TLS int Poly_Impl::remove_space_dims_percentage = 20;

} // namespace pplite
