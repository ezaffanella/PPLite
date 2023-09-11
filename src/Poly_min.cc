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

#include "pplite-config.h"

#include "Poly.hh"
#include "Poly_templ.hh"
#include "Poly_min.hh"
#include "Var.hh"

#include <algorithm>
#include <functional>
#include <utility>

namespace pplite {

namespace detail {

void
set_universe(const dim_type space_dim,
             Poly_Impl::Sys<Cons>& cs,
             Poly_Impl::Sys<Gens>& gs,
             Sat& sat_c) {
  assert(cs.empty() && gs.empty());
  assert(sat_c.num_rows() == 0 && sat_c.num_cols() == 0);
  // Add the strict positivity constraint.
  cs.sk_rows.push_back(Con::zero_dim_positivity());
  // Add space_dim lines.
  if (space_dim > 0) {
    auto& lines = gs.sing_rows;
    lines.reserve(space_dim);
    for (dim_type i = 0; i < space_dim; ++i)
      lines.push_back(line(Var(i)));
  }
  // Add the point at the origin of the vector space.
  gs.sk_rows.push_back(point());
  // The point does not saturate the positivity constraint.
  sat_c.resize(1, 1);
  sat_c[0].set(0);
}

void
init_dd(const dim_type space_dim,
        Poly_Impl::Sys<Gens>& src,
        Poly_Impl::Sys<Gens>& pending,
        Poly_Impl::Sys<Cons>& dst,
        Sat& sat_g) {
  assert(src.empty() && !pending.empty() && dst.empty());
  assert(sat_g.num_rows() == 0 && sat_g.num_cols() == 0);
  // Search for a point.
  auto iter_point
    = std::find_if(pending.sk_rows.begin(), pending.sk_rows.end(),
                   std::mem_fn(&Gen::is_point));
  if (iter_point != pending.sk_rows.end()) {
    // Found a point: this will be moved into src.
    // Before moving it, remap the pending supports.
    auto index = std::distance(pending.sk_rows.begin(), iter_point);
    Index_Set tbr(index);
    for (auto& ns : pending.ns_rows)
      ns.remove_all(tbr);
    // Now move the point found into src.
    src.sk_rows.push_back(std::move(*iter_point));
    pending.sk_rows.erase(iter_point);
  } else {
    // No point found: materialize a non-skel point.
    assert(!pending.ns_rows.empty());
    // Choosing last one (simpler to remove later).
    const auto& ns = pending.ns_rows.back();
    Gen mat_point = materialize(ns, pending.sk_rows);
    src.sk_rows.push_back(std::move(mat_point));
    pending.ns_rows.pop_back();
  }
  // Add space_dim equalities.
  if (space_dim > 0) {
    const auto& p = src.sk_rows[0];
    const auto& p_div = p.divisor();
    auto& equals = dst.sing_rows;
    equals.reserve(space_dim);
    for (dim_type i = 0; i < space_dim; ++i) {
      Linear_Expr expr;
      expr.set_space_dim(i+1);
      expr.set(i, p_div);
      Con eq(std::move(expr), -p.coeff(Var(i)), Con::EQUALITY);
      equals.push_back(std::move(eq));
    }
  }
  // Add the positivity constraint.
  dst.sk_rows.push_back(Con::zero_dim_positivity());
  // The point does not saturate the positivity constraint.
  sat_g.resize(1, 1);
  sat_g[0].set(0);
}

} // namespace detail

// Initialize static (thread local) data member.
PPLITE_TLS dim_type Poly_Impl::minimize_filter_threshold = not_a_dim();

} // namespace pplite

