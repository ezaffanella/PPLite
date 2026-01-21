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

#include "globals.hh"
#include "support_utils.hh"
#include "utils.hh"
#include "Bits.hh"
#include "Con.hh"
#include "Gen.hh"
#include "Sat.hh"
#include "Poly_min.hh"
#include "Var.hh"

#include <cassert>
#include <utility>

namespace pplite {
namespace detail {

template <typename Sys>
inline void
detect_implicit_singular(Sys& sys, Sat& sat) {
  auto& sg_rows = sys.sing_rows;
  auto& sk_rows = sys.sk_rows;
  auto& ns_rows = sys.ns_rows;
  assert(num_rows(sk_rows) == sat.num_rows());
  // Identify and move implicit singular rows.
  Dims indices;
  for (auto i : range(sat.num_rows())) {
    if (sat[i].empty()) {
      indices.push_back(i);
      sg_rows.push_back(std::move(sk_rows[i]));
      make_singular(sg_rows.back());
    }
  }
  // Delete them from sk_rows and sat (keeping ns_rows consistent).
  remove_rows(indices, sk_rows, ns_rows, sat);
}

template <typename Sys>
inline void
simplify_singular(Sys& sys, Sat& sat) {
  detect_implicit_singular(sys, sat);
  // Remove redundant singular elements.
  gauss(sys.sing_rows);
}

Bits
sat_positivity(const Gens& gs) {
  Bits row;
  for (auto i : bwd_index_range(gs)) {
    if (gs[i].is_point_or_closure_point())
      row.set(i);
  }
  return row;
}

void
normalize_strict_positivity(Cons& sk_src, NS_Rows& ns_src,
                            const Gens& sk_dst, Sat& sat_dst) {
  // Compute sat row for the positivity constraint.
  Bits sat_pos = sat_positivity(sk_dst);

  auto is_strict_positivity = [&sk_src, &sat_pos, &sat_dst](dim_type i)
    {
      return sk_src[i].is_strict_inequality()
      && sat_pos == sat_dst[i];
    };

  // Find the first skel strict ineq that saturates no closure point.
  const auto last = num_rows(sk_src);
  auto first = last;
  for (auto i : range(last)) {
    if (is_strict_positivity(i)) {
      first = i;
      break;
    }
  }
  if (first == last)
    return;

  // Replace it with a copy of the (strict) positivity constraint.
  sk_src[first] = Con::zero_dim_positivity();

  // Other occurrences (if any) are removed.
  Dims tbr;
  for (auto i : range(first+1, last)) {
    if (is_strict_positivity(i))
      tbr.push_back(i);
  }
  remove_rows(tbr, sk_src, ns_src, sat_dst);
}

inline void
normalize_strict_positivity(Gens&, NS_Rows&, const Cons&, Sat&) {
  // Nothing to be done.
}

template <typename SK_Rows>
inline void
simplify_skel_quick_check(const dim_type min_sat,
                          SK_Rows& sk_rows, NS_Rows& ns_rows, Sat& sat) {
  const auto num_cols = sat.num_cols();
  Dims redundant;
  for (auto i : index_range(sk_rows)) {
    // FIXME: quando si lavora sui generatori, scartare a priori tutti
    // i punti porta a inefficienze; dovremo farlo solo per gli NNC.
    if (is_strict_ineq_or_point(sk_rows[i].type()))
      continue;
    const auto num_sat = num_cols - sat[i].count_ones();
    assert(num_sat < num_cols);
    if (num_sat < min_sat)
      redundant.push_back(i);
  }
  remove_rows(redundant, sk_rows, ns_rows, sat);
}

template <typename Row_Type>
inline bool
is_made_redundant_by(Row_Type type_i, const Bits& sat_i,
                     Row_Type type_j, const Bits& sat_j) {
  assert(!is_singular(type_i));
  assert(!is_singular(type_j));
  if (type_i == type_j)
    return subset_eq(sat_j, sat_i);
  // Note: this check is ad-hoc for NNC case.
  if (is_nonstrict_ineq_or_closure_point(type_i)
      && is_strict_ineq_or_point(type_j)) {
    return subset_eq(sat_j, sat_i);
  }
  // i is not made redundant by j.
  return false;
}

template <typename SK_Rows>
inline void
simplify_skel_full_check(SK_Rows& sk_rows, NS_Rows& ns_rows, Sat& sat) {
  Index_Set redundant;
  for (auto i : index_range(sk_rows)) {
    const auto type_i = sk_rows[i].type();
    const auto& sat_i = sat[i];
    for (auto j : index_range(sk_rows)) {
      if (i == j || (i > j && redundant.test(j)))
        continue;
      const auto type_j = sk_rows[j].type();
      const auto& sat_j = sat[j];
      if (is_made_redundant_by(type_i, sat_i, type_j, sat_j)) {
        redundant.set(i);
        break;
      }
    }
  }
  remove_rows(redundant, sk_rows, ns_rows, sat);
}

template <typename Row_Type>
inline bool
is_made_nonskel_by(Row_Type type_i, const Bits& sat_i,
                   Row_Type type_j, const Bits& sat_j) {
  assert(!is_singular(type_i));
  assert(!is_singular(type_j));
  return is_strict_ineq_or_point(type_i)
    && is_nonstrict_ineq_or_closure_point(type_j)
    && subset_ne(sat_j, sat_i);
}

template <typename SK_Rows>
inline void
simplify_non_skel(SK_Rows& sk_rows, NS_Rows& ns_rows,
                  Sat& sat, Sat& sat_tr) {
  Dims non_skel;
  for (auto i : index_range(sk_rows)) {
    const auto type_i = sk_rows[i].type();
    const auto& sat_i = sat[i];
    for (auto j : index_range(sk_rows)) {
      if (i == j)
        continue;
      const auto type_j = sk_rows[j].type();
      const auto& sat_j = sat[j];
      if (is_made_nonskel_by(type_i, sat_i, type_j, sat_j)) {
        non_skel.push_back(i);
        break;
      }
    }
  }

  if (non_skel.size() > 0) {
    // Move non-skel elements into ns_rows.
    Sat sat_ns(non_skel.size(), sat.num_cols());
    for (auto i : index_range(non_skel)) {
      auto ns_index = non_skel[i];
      assert(is_strict_ineq_or_point(sk_rows[ns_index].type()));
      using std::swap;
      swap(sat_ns[i], sat[ns_index]);
    }
    remove_rows(non_skel, sk_rows, ns_rows, sat);

    sat_tr = sat.transpose();
    const auto sat_tr_cols = sat_tr.num_cols();
    for (auto& row : sat_ns.impl().rows) {
      // This is one half of support_closure().
      row.complement_until(sat.num_cols());
      Index_Set iset(std::move(row));
      row = sat_all(iset, sat_tr);
      row.complement_until(sat_tr_cols);
      ns_rows.push_back(Index_Set(std::move(row)));
    }
  }

  // Detect non-minimal support sets.
  Index_Set redundant;
  for (auto i : index_range(ns_rows)) {
    if (redundant.test(i))
      continue;
    for (auto j : index_range(ns_rows)) {
      if (i == j || redundant.test(j))
        continue;
      if (subset_eq(ns_rows[i], ns_rows[j]))
        redundant.set(j);
    }
  }

  // Detect support sets containing a skel element that could
  // be a non-skel element (i.e., point or strict inequality).
  Index_Set skel_nonskel;
  for (auto i : index_range(sk_rows)) {
    if (is_strict_ineq_or_point(sk_rows[i].type()))
      skel_nonskel.set(i);
  }
  for (auto i : index_range(ns_rows)) {
    if (redundant.test(i))
      continue;
    if (ns_rows[i].intersects(skel_nonskel))
      redundant.set(i);
  }
  erase_using_sorted_indices(ns_rows, redundant);

  // Check for ns supported by a single skel element.
  promote_singletons(ns_rows, sk_rows);
}

////////////////////////////////////////////////////////////////////////

// Template definition.
template <typename Src_Sys, typename Dst_Sys>
void
simplify(const dim_type space_dim,
         Src_Sys& src, Sat& sat_src, Sat& sat_dst,
         const Dst_Sys& dst) {
  assert(sat_src.num_rows() == num_rows(dst.sk_rows));
  assert(sat_src.num_cols() == num_rows(src.sk_rows));

  sat_dst = sat_src.transpose();
  simplify_singular(src, sat_dst);
  // Maybe restore sat_src (possibly made invalid by previous call).
  if (sat_dst.num_rows() != sat_src.num_cols())
    sat_src = sat_dst.transpose();

  auto& sk_src = src.sk_rows;
  auto& ns_src = src.ns_rows;

  // Recompute support closure, so as to replace the redundant
  // supporting elements by those making them redundant.
  for (auto& ns : ns_src)
    support_closure(ns, sat_dst, sat_src);

  // Next call is specific for Src_Sys = Con_Sys
  // (it does nothing at all when Src_Sys = Gen_Sys).
  normalize_strict_positivity(sk_src, ns_src, dst.sk_rows, sat_dst);

  const auto num_src_sing = num_rows(src.sing_rows);
  assert(num_src_sing <= space_dim);
  const auto num_dst_sing = num_rows(dst.sing_rows);
  assert(num_dst_sing <= space_dim);
  const auto min_sat = space_dim - num_src_sing - num_dst_sing;

  simplify_skel_quick_check(min_sat, sk_src, ns_src, sat_dst);
  simplify_skel_full_check(sk_src, ns_src, sat_dst);
  simplify_non_skel(sk_src, ns_src, sat_dst, sat_src);
  assert(check_redundant_ns_by_ns(ns_src));
  // Maybe restore sat_src (possibly made invalid by previous calls).
  if (sat_dst.num_rows() != sat_src.num_cols())
    sat_src = sat_dst.transpose();

  // Redundancies removed: now perform back substitution.
  back_substitute(src.sing_rows, src.sk_rows);

  maybe_ensure_strict_pos(sk_src, ns_src);
  assert(check_redundant_ns_by_ns(ns_src));
  assert(maybe_check_efc(ns_src, sk_src, dst.sk_rows, sat_dst, sat_src));
}

// Explicit instantiations.
template
void
simplify<Con_Sys, Gen_Sys>(dim_type, Con_Sys&, Sat&, Sat&, const Gen_Sys&);
template
void
simplify<Gen_Sys, Con_Sys>(dim_type, Gen_Sys&, Sat&, Sat&, const Con_Sys&);

} // namespace detail

} // namespace pplite
