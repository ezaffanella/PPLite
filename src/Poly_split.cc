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

#include "globals.hh"
#include "utils.hh"
#include "support_utils.hh"
#include "Bits.hh"
#include "Con.hh"
#include "Gen.hh"
#include "Index_Partition.hh"
#include "Integer.hh"
#include "Integer.hh"
#include "Linear_Expr.hh"
#include "Poly.hh"
#include "Poly_templ.hh"
#include "Poly_min.hh"
#include "Poly_Rel.hh"
#include "Sat.hh"
#include "Scalar_Prod.hh"

#include <algorithm>
#include <cassert>
#include <utility>
#include <vector>

namespace pplite {
namespace detail {

namespace {

using Split_Result = std::pair<Conv_Iter_Result, Conv_Iter_Result>;

inline bool
has_no_eq_point(const Ranges& sk_ranges,
                const Gens& sk, const NS_Rows& ns_rows) {
  Index_Set sk_eq = sk_ranges_to_eq_index_set(sk_ranges);
  for (const auto& g : sk_eq)
    if (sk[g].is_point())
      return false;
  return !any_of(ns_rows, [&sk_eq](const Index_Set& ns) {
                            return subset_eq(ns, sk_eq);
                          });
}

inline bool
has_no_eq_point(const Ranges& sk_ranges, const Gens& sk) {
  NS_Rows dummy;
  return has_no_eq_point(sk_ranges, sk, dummy);
}

/*******************************************************/
// Process skeleton

// Combines the skeleton generator in SK+ with those in SK-.
// Creates new elements in sk_strict and, if nnc_split, also in sk_nonstrict.
// Updates only the `sat_strict': its copy is postponed.
void
Q_split_combine(Ranges& sk_ranges, Gens& sk_strict, Gens& sk_nonstrict,
                Integers& sp, Sat& sat_strict,
                const bool nnc_split,
                const dim_type min_sat, const dim_type max_sat) {
  const dim_type old_num_sk_rows = sk_strict.size();

  static PPLITE_TLS Dims sk_pos, sk_neg;
  static PPLITE_TLS Dims sat_ones;
  sk_pos.clear();
  sk_neg.clear();
  sat_ones.clear();
  for (auto i : index_range(sk_ranges)) {
    if (sk_ranges[i] == Range::EQ) {
      sat_ones.push_back(0);
      continue;
    }
    sat_ones.push_back(sat_strict[i].count_ones());
    if (sk_ranges[i] == Range::POS)
      sk_pos.push_back(i);
    else
      sk_neg.push_back(i);
  }

  for (auto pos : sk_pos) {
    for (auto neg : sk_neg) {
      // Exploit saturation information to perform adjacency tests.
      // Note: delay creation of `new_row' and `new_satrow'.
      const dim_type new_satrow_ones
        = sat_strict[pos].count_ones_in_union(sat_strict[neg]);
      if (quick_non_adj_test(min_sat, max_sat, new_satrow_ones))
        continue;
      Bits new_satrow = Bits::get_union(sat_strict[pos], sat_strict[neg]);
      if (not quick_adj_test(pos, neg, new_satrow_ones, sat_ones) &&
          not combinatorial_adj_test(0, old_num_sk_rows,
                                     pos, neg, new_satrow, sat_strict))
        continue;

      // Here `pos' and `neg' are adjacent: update `sat_strict'.
      sat_strict.add_row(std::move(new_satrow));

      // Add the new row to `sk_strict'
      auto new_row_strict = sk_strict[neg];
      // Combination can be done once.
      combine_row(new_row_strict, sk_strict[pos], sp[neg], sp[pos]);
      new_row_strict.set_type(combine_type(sk_strict[pos].type(),
                                           sk_strict[neg].type(),
                                           nnc_split));
      if (nnc_split) {
        // Copy the new element
        auto new_row_nonstrict = new_row_strict;
        // Set the right type, according to the side.
        new_row_nonstrict.set_type(combine_type(sk_strict[neg].type(),
                                                sk_strict[pos].type(),
                                                false));
        sk_nonstrict.push_back(std::move(new_row_nonstrict));
      }
      sk_strict.push_back(std::move(new_row_strict));

      assert(num_rows(sk_strict) == sat_strict.num_rows());
    } // end loop on sk_neg
  } // end loop on sk_pos

  sp.resize(sk_strict.size());
  sk_ranges.resize(sk_strict.size(), Range::EQ);
}

/**
 * Computes the split sk systems, copying the original one.
 * NB: sk_[non]strict maskes sense iff [non]strict result is not EMPTY.
 * If a component is EMPTY, it will be the caller to set the empty DD.
 */
Split_Result
Q_split_skel(dim_type space_dim, dim_type num_lines, dim_type num_eq,
             Gens& sk_strict, Gens& sk_nonstrict,
             Ranges& sk_ranges, const Con& beta_strict,
             const NS_Rows& ns_strict, Sat& sat_strict) {
  assert(beta_strict.is_inequality());
  const bool nnc_split = beta_strict.is_strict_inequality();
  bool empty_strict = false;
  bool empty_nonstrict = false;

  Integers scal_prods;
  sp::compute_prods(beta_strict, sk_strict, scal_prods);

  // Create partition for skeleton rows.
  // The partition is relative to the strict side.
  sk_partition(scal_prods, sk_ranges);
  Index_Partition<Range::NEG> sk_neg(sk_ranges);
  Index_Partition<Range::EQ> sk_eq(sk_ranges);
  Index_Partition<Range::POS> sk_pos(sk_ranges);

  if (!sk_neg.empty() && !sk_pos.empty()) {
    // Here SK+ and SK- are non-empty.
    // Combine SK+ and SK- *adjacent* generators.
    const dim_type min_sat = space_dim - num_lines - 1;
    const dim_type max_sat = sat_strict.num_cols() + num_eq;
    if (nnc_split) {
      // Copy sk_strict into sk_nonstrict
      sk_nonstrict = sk_strict;
    }
    Q_split_combine(sk_ranges, sk_strict, sk_nonstrict,
                    scal_prods, sat_strict,
                    nnc_split, min_sat, max_sat);
  } else {
    // Here SK+ or SK- is empty.
    const bool no_eq_point = has_no_eq_point(sk_ranges, sk_strict, ns_strict);
    if (sk_neg.empty()) {
      if (sk_eq.empty() || no_eq_point)
        empty_nonstrict = true;
      if (nnc_split && sk_pos.empty()) {
        // strict side is only in eq, but it is strict.
        empty_strict = true;
      }
      // The non-skel handling is postponed.
    }
    else {
      assert(sk_pos.empty());
      if (nnc_split || sk_eq.empty() || no_eq_point)
        empty_strict = true;
    }

    // Copy original gs system only if both sides are non-empty.
    // If !nnc_split the copy is postponed, since there could be
    // shared promotions of the same ns elements.
    if (!empty_strict && !empty_nonstrict && nnc_split)
      sk_nonstrict = sk_strict;
  }

  assert(not (empty_strict && empty_nonstrict));
  if (empty_strict)
    return { Conv_Iter_Result::EMPTY, Conv_Iter_Result::REDUNDANT };
  if (empty_nonstrict)
    return { Conv_Iter_Result::REDUNDANT, Conv_Iter_Result::EMPTY };

  // Update sat_strict cols.
  sat_strict.add_cols(1);
  update_last_sat_column(sat_strict, sk_pos);
  update_last_sat_column(sat_strict, sk_neg);
  return { Conv_Iter_Result::OK, Conv_Iter_Result::OK };
}

/*******************************************************/
// Process non-skeleton

/* Move-ns method for NNC split */
// Computes the closures once and differentiates only the projection phase.
// Checks the supports to be promoted only in the nonstrict side.
void
split_move_ns(NS_Rows& ns_strict, NS_Rows& ns_nonstrict,
              Ranges& ns_ranges_strict, Ranges& ns_ranges_nonstrict,
              Gens& sk_nonstrict,
              const Index_Set& iset_neg,
              const Index_Set& iset_eq,
              const Index_Set& iset_pos,
              NS_Rows& nsstar_strict, NS_Rows& nsstar_nonstrict,
              const Sat& sat_c, Sat& sat_g) {
  // Originally ns_ranges are the same.
  Dims to_be_del;
  Index_Partition<Range::POS_NEG> ns_pos_neg(ns_ranges_strict);
  for (auto i : ns_pos_neg) {
    assert(ns_ranges_strict[i] == Range::POS_NEG);
    Index_Set& supp_strict = ns_strict[i];
    // Checkme: closures can be done once, since sat info are shared.
    support_closure(supp_strict, sat_c, sat_g);

    Index_Set supp_nonstrict = supp_strict;
    // Differentiate only the projection phase.
    supp_projection(supp_strict, true, iset_neg, iset_eq);
    supp_projection(supp_nonstrict, false, iset_pos, iset_eq);
    // No need to check if it has to be promoted
    assert(supp_strict.size() > 1);
    add_to_minimal_set_of_ns(supp_strict, nsstar_strict);

    // Check if it can be moved into SK.
    if (supp_nonstrict.size() == 1) {
      // Change its only supporter's type into point.
      const auto supporter_index = *(supp_nonstrict.begin());
      promote_singleton(sk_nonstrict[supporter_index]);
      to_be_del.push_back(i);
      // Consider next ns in pos_neg
      continue;
    }
    assert(supp_nonstrict.size() > 1);
    add_to_minimal_set_of_ns(std::move(supp_nonstrict), nsstar_nonstrict);
  }
  // Remove ns that has been moved into SK.
  erase_using_sorted_indices(ns_nonstrict, to_be_del);
  erase_using_sorted_indices(ns_ranges_nonstrict, to_be_del);
}

/** Create_ns method for NNC split. **/
// Add to NS*strict:
// - ns* from NEG projected in POS,
// - ns* from EQ projected in POS.
// Add to NS*nonstrict:
// - ns* from NEG projected in EQ,
// - ns* from POS projected in EQ.
void
split_create_ns(NS_Rows& ns_strict, NS_Rows& ns_nonstrict,
                Ranges& ns_ranges_strict, Ranges& ns_ranges_nonstrict,
                Gens& sk_strict, Gens& sk_nonstrict, const Ranges& sk_ranges,
                const Index_Set& iset_pos,
                const Index_Set& iset_eq,
                const Index_Set& iset_neg,
                NS_Rows& nsstar_strict, NS_Rows& nsstar_nonstrict,
                const Sat& sat_c, const Sat& sat_g,
                dim_type min_sat, dim_type max_sat,
                bool should_create_nonstrict) {
  // enumerate_faces from NEG can be done once:
  // Store the computed set, reproject it in EQ if needed.
  NS_Rows set_from_neg;
  create_ns_from_Q_with_face_enum<Range::NEG>(Range::POS,
                                              sk_strict, ns_strict,
                                              sk_ranges, ns_ranges_strict,
                                              set_from_neg,
                                              iset_pos, iset_eq, iset_neg,
                                              sat_c, sat_g,
                                              true);
  // Here `set_from_neg' stores the created ns projected in POS.
  for (auto& ns_star : set_from_neg) {
    // Add it to nsstar_strict (it is copied).
    add_to_minimal_set_of_ns(ns_star, nsstar_strict);
    if (should_create_nonstrict) {
      // Reproject it in EQ
      supp_projection(ns_star, false, iset_neg, iset_eq);
      // Add it to nsstar_nonstrict too.
      add_to_minimal_set_of_ns(std::move(ns_star), nsstar_nonstrict);
    }
  }

  // Add to nsstar_strict new ns from EQ.
  create_ns_from_eq_points(sk_strict, ns_strict, sk_ranges, ns_ranges_strict,
                           nsstar_strict, iset_pos, iset_eq, iset_neg,
                           sat_c, sat_g, min_sat, max_sat);

  if (should_create_nonstrict) {
    // Add to nsstar_nonstrict new ns from POS.
    create_ns_from_Q_with_face_enum<Range::POS>(Range::NEG,
                                                sk_nonstrict, ns_nonstrict,
                                                sk_ranges, ns_ranges_nonstrict,
                                                nsstar_nonstrict,
                                                iset_pos, iset_eq, iset_neg,
                                                sat_c, sat_g,
                                                false);
  }
}

void
split_non_skel(NS_Rows& ns_strict, NS_Rows& ns_nonstrict,
               Gens& sk_strict, Gens& sk_nonstrict,
               const Ranges& sk_ranges,
               const Sat& sat_c, Sat& sat_g,
               dim_type space_dim, dim_type num_lines, dim_type num_eq) {
  NS_Rows nsstar_strict, nsstar_nonstrict;
  Ranges ns_ranges_strict, ns_ranges_nonstrict;
  // Create partition for non-skeleton rows.
  Index_Set iset_neg, iset_eq, iset_pos;
  sk_ranges_to_index_sets(sk_ranges, iset_neg, iset_eq, iset_pos);
  ns_partition(ns_strict, iset_neg, iset_pos, ns_ranges_strict);

  // Copy partition, since they will follow different ns elements.
  // NB: the real nonstrict partition should be its opposite,
  // it is not complemented for efficiency.
  ns_nonstrict = ns_strict;
  ns_ranges_nonstrict = ns_ranges_strict;
  nsstar_strict.clear();
  nsstar_nonstrict.clear();

  bool sat_is_up_to_date = false;
  Index_Partition<Range::POS_NEG> ns_pos_neg(ns_ranges_strict);
  if (!ns_pos_neg.empty()) {
    sat_g = sat_c.transpose();
    sat_is_up_to_date = true;
    split_move_ns(ns_strict, ns_nonstrict,
                  ns_ranges_strict, ns_ranges_nonstrict,
                  sk_nonstrict,
                  iset_neg, iset_eq, iset_pos,
                  nsstar_strict, nsstar_nonstrict,
                  sat_c, sat_g);
  }

  // Check if we need to call `create_ns' for the non-strict side too.
  const bool should_create_nonstrict =
    can_have_nonskel(iset_eq, sk_nonstrict);
  if (!sat_is_up_to_date) {
    sat_g = sat_c.transpose();
    sat_is_up_to_date = true;
  }
  const dim_type min_sat = space_dim - num_eq - 1;
  const dim_type max_sat = sat_c.num_cols() + num_lines;
  split_create_ns(ns_strict, ns_nonstrict,
                  ns_ranges_strict, ns_ranges_nonstrict,
                  sk_strict, sk_nonstrict, sk_ranges,
                  iset_pos, iset_eq, iset_neg,
                  nsstar_strict, nsstar_nonstrict,
                  sat_c, sat_g, min_sat, max_sat,
                  should_create_nonstrict);
  // Add new ns.
  if (sat_is_up_to_date) {
    new_ns_from_set(sk_strict, ns_strict, nsstar_strict,
                    sk_ranges, ns_ranges_strict, true);
    new_ns_from_set(sk_nonstrict, ns_nonstrict, nsstar_nonstrict,
                    sk_ranges, ns_ranges_nonstrict, false);
  }

  // Final update
  std::replace(ns_ranges_strict.begin(), ns_ranges_strict.end(),
               Range::POS_NEG, Range::NEG);
  std::replace(ns_ranges_strict.begin(), ns_ranges_strict.end(),
               Range::EQ, Range::NEG);
  std::replace(ns_ranges_nonstrict.begin(), ns_ranges_nonstrict.end(),
               Range::POS_NEG, Range::POS);
  Index_Partition<Range::NEG> ns_neg(ns_ranges_strict);
  Index_Partition<Range::POS> ns_pos(ns_ranges_nonstrict);
  erase_using_sorted_indices(ns_strict, ns_neg);
  erase_using_sorted_indices(ns_nonstrict, ns_pos);
}

/*******************************************************/
// Non-skel handling in CLOSED split.

/* Overloaded move-ns method for CLOSED split */
// Compute the closures and project in EQ.
// Check the supports to be promoted.
void
split_move_ns(NS_Rows& ns_rows, Ranges& ns_ranges, Gens& sk_rows,
              const Index_Set& iset_neg,
              const Index_Set& iset_eq,
              NS_Rows& nsstar,
              const Sat& sat_c, Sat& sat_g) {
  // Originally ns_ranges are the same.
  Dims to_be_del;
  Index_Partition<Range::POS_NEG> ns_pos_neg(ns_ranges);
  for (auto i : ns_pos_neg) {
    assert(ns_ranges[i] == Range::POS_NEG);
    Index_Set& supp = ns_rows[i];
    support_closure(supp, sat_c, sat_g);
    supp_projection(supp, false, iset_neg, iset_eq);

    // Check if it can be moved into SK.
    if (supp.size() == 1) {
      // Change its only supporter's type into point.
      const auto supporter_index = *(supp.begin());
      promote_singleton(sk_rows[supporter_index]);
      to_be_del.push_back(i);
      // Consider next ns in pos_neg
      continue;
    }
    add_to_minimal_set_of_ns(supp, nsstar);
  }
  // Remove ns that has been moved into SK.
  erase_using_sorted_indices(ns_rows, to_be_del);
  erase_using_sorted_indices(ns_ranges, to_be_del);
}

/* Overloaded create_ns method for CLOSED split. */
// Add to NS*:
// - ns* from NEG projected in EQ,
// - ns* from POS projected in EQ.
void
split_create_ns(NS_Rows& ns_rows, Ranges& ns_ranges,
                Gens& sk_rows, const Ranges& sk_ranges,
                const Index_Set& iset_pos,
                const Index_Set& iset_eq,
                const Index_Set& iset_neg,
                NS_Rows& nsstar,
                const Sat& sat_c, const Sat& sat_g) {
  create_ns_from_Q_with_face_enum<Range::NEG>(Range::POS,
                                              sk_rows, ns_rows,
                                              sk_ranges, ns_ranges,
                                              nsstar,
                                              iset_pos, iset_eq, iset_neg,
                                              sat_c, sat_g,
                                              false);
  create_ns_from_Q_with_face_enum<Range::POS>(Range::NEG,
                                              sk_rows, ns_rows,
                                              sk_ranges, ns_ranges,
                                              nsstar,
                                              iset_pos, iset_eq, iset_neg,
                                              sat_c, sat_g,
                                              false);
}

/* Overloaded version for CLOSED split */
// Computes only one NS* and copies the whole NS sys only in the end.
void
split_non_skel(NS_Rows& ns_strict, NS_Rows& ns_nonstrict,
               Gens& sk_rows, const Ranges& sk_ranges,
               const Sat& sat_c, Sat& sat_g) {
  NS_Rows nsstar_strict;
  Ranges ns_ranges_strict;
  // Create partition for non-skeleton rows.
  Index_Set iset_neg, iset_eq, iset_pos;
  sk_ranges_to_index_sets(sk_ranges, iset_neg, iset_eq, iset_pos);
  ns_partition(ns_strict, iset_neg, iset_pos, ns_ranges_strict);
  nsstar_strict.clear();

  bool sat_is_up_to_date = false;
  Index_Partition<Range::POS_NEG> ns_pos_neg(ns_ranges_strict);
  if (!ns_pos_neg.empty()) {
    sat_g = sat_c.transpose();
    sat_is_up_to_date = true;
    split_move_ns(ns_strict, ns_ranges_strict, sk_rows,
                  iset_neg, iset_eq,
                  nsstar_strict, sat_c, sat_g);
  }

  // Check if we need to call `create_ns'.
  const bool should_create_nonstrict =
    can_have_nonskel(iset_eq, sk_rows);
  if (should_create_nonstrict) {
    if (!sat_is_up_to_date) {
      sat_g = sat_c.transpose();
      sat_is_up_to_date = true;
    }
    split_create_ns(ns_strict, ns_ranges_strict,
                    sk_rows, sk_ranges,
                    iset_pos, iset_eq, iset_neg,
                    nsstar_strict,
                    sat_c, sat_g);
  }
  // Add new ns.
  if (sat_is_up_to_date) {
    new_ns_from_set(sk_rows, ns_strict, nsstar_strict,
                    sk_ranges, ns_ranges_strict, true);
  }

  // Copy NS systems.
  Ranges ns_ranges_nonstrict = ns_ranges_strict;
  ns_nonstrict = ns_strict;

  // Final update
  std::replace(ns_ranges_strict.begin(), ns_ranges_strict.end(),
               Range::POS_NEG, Range::NEG);
  std::replace(ns_ranges_nonstrict.begin(), ns_ranges_nonstrict.end(),
               Range::POS_NEG, Range::POS);
  Index_Partition<Range::NEG> ns_neg(ns_ranges_strict);
  Index_Partition<Range::POS> ns_pos(ns_ranges_nonstrict);
  erase_using_sorted_indices(ns_strict, ns_neg);
  erase_using_sorted_indices(ns_nonstrict, ns_pos);
}

/*******************************************************/
// Process violating line

void
split_violating_singular(dim_type ex,
                         Gen_Sys& sys_strict,
                         Gen_Sys& sys_nonstrict,
                         Sat& sat_c_strict, Sat& sat_c_nonstrict,
                         const Con& beta_strict) {
  auto& sg_strict = sys_strict.sing_rows;
  auto& sk_strict = sys_strict.sk_rows;
  auto& ns_strict = sys_strict.ns_rows;
  auto& sk_nonstrict = sys_nonstrict.sk_rows;

  const bool nnc_split = beta_strict.is_strict_inequality();

  // Map the ex_strict in a ray
  auto& ex_strict = sg_strict[ex];
  Integer sp_ex;
  sp::assign(sp_ex, beta_strict, ex_strict);
  make_non_singular(ex_strict, sp_ex);
  // Map all other geometric rows into SK=: can be done once.
  combine_with_ex_sing(ex, sp_ex, beta_strict, sg_strict, sk_strict);
  // Move the ex singular to sk
  const auto new_ex = num_rows(sk_strict);
  sk_strict.push_back(std::move(ex_strict));
  sg_strict.erase(sg_strict.begin() + ex);
  // Copy sys: change ray orientation for the non-strict side.
  sys_nonstrict = sys_strict;
  auto& ex_nonstrict = sk_nonstrict[num_rows(sk_nonstrict) - 1];
  assert(ex_nonstrict == sk_strict[new_ex]);
  neg_assign(ex_nonstrict.linear_expr());

  // Update sat_c : can be done once.
  assert(beta_strict.is_inequality());
  sat_c_strict.add_cols(1);
  Bits sat_ex;
  sat_ex.set(sat_c_strict.num_cols() - 1);
  sat_c_strict.add_row(std::move(sat_ex));
  // Copy sat: sat info are shared.
  sat_c_nonstrict = sat_c_strict;

  // Here sk_nonstrict and ns_nonstrict are done.

  if (nnc_split) {
    // Handle ns_strict and sk_strict eq points.
    Ranges sk_ranges, ns_ranges;
    sk_ranges.assign(sk_strict.size(), Range::EQ);
    sk_ranges[new_ex] = Range::POS;
    ns_ranges.assign(ns_strict.size(), Range::EQ);
    create_ns_with_ex_sing(sk_strict, ns_strict,
                           sk_ranges, ns_ranges, new_ex);
  }
}

/********************************************************/

/**
 * Split with RATIONAL VARIABLES.
 * This is the usual conversion algorithm, only
 * updating gen_sys_then and gen_sys_else differently.
 */
Split_Result
Q_split_gen_sys(dim_type space_dim, dim_type num_eq,
                Gen_Sys& gen_sys_strict, Gen_Sys& gen_sys_nonstrict,
                Sat& sat_c_strict, Sat& sat_c_nonstrict,
                Sat& sat_g_strict,
                const Con& beta_strict, bool nnc_poly) {
  // If `nnc_split' then the "strict" side is actually strict.
  const bool nnc_split = beta_strict.is_strict_inequality();
  assert(!nnc_split || nnc_poly);

  Gens& sg_strict = gen_sys_strict.sing_rows;
  Gens& sk_strict = gen_sys_strict.sk_rows;
  Gens& sk_nonstrict = gen_sys_nonstrict.sk_rows;
  NS_Rows& ns_strict = gen_sys_strict.ns_rows;
  NS_Rows& ns_nonstrict = gen_sys_nonstrict.ns_rows;

  // Check if all singular rows saturate strict_beta.
  const dim_type ex = check_sing_rows(beta_strict, sg_strict);
  if (ex < num_rows(sg_strict)) {
    split_violating_singular(ex, gen_sys_strict,
                             gen_sys_nonstrict,
                             sat_c_strict, sat_c_nonstrict,
                             beta_strict);
    return { Conv_Iter_Result::OK, Conv_Iter_Result::OK };
  }

  // Here all lines saturate the constraint.
  // Copy sing sys
  gen_sys_nonstrict.sing_rows = sg_strict;
  dim_type num_lines = num_rows(sg_strict);

  // Split the skeleton part.
  // Copies the sk-sys only if nnc_split and if both sides are non-empty.
  // Updates only `sat_strict': the copy is postponed.
  Ranges sk_ranges;
  auto [res_strict, res_nonstrict]
    = Q_split_skel(space_dim, num_lines, num_eq,
                   sk_strict, sk_nonstrict,
                   sk_ranges, beta_strict,
                   ns_strict, sat_c_strict);

  if (res_strict == Conv_Iter_Result::EMPTY ||
      res_nonstrict == Conv_Iter_Result::EMPTY)
    return { res_strict, res_nonstrict };

  // Non-skel handling is required only in NNC polyhedra.
  if (nnc_poly) {
    if (nnc_split) {
      split_non_skel(ns_strict, ns_nonstrict,
                     sk_strict, sk_nonstrict, sk_ranges,
                     sat_c_strict, sat_g_strict,
                     space_dim, num_lines, num_eq);
    } else {
      split_non_skel(ns_strict, ns_nonstrict, sk_strict, sk_ranges,
                     sat_c_strict, sat_g_strict);
    }
  }

  // Copy sat: all saturation info is shared.
  sat_c_nonstrict = sat_c_strict;
  // Note: we do NOT copy sat_g_strict into sat_g_nonstrict,
  // because these will be anyway rewritten by calls to simplify.
  if (not nnc_split)
    // Copy SK sys since it has been computed only once.
    sk_nonstrict = sk_strict;

  // Final SK update.
  Index_Partition<Range::NEG> sk_neg(sk_ranges);
  Index_Partition<Range::EQ> sk_eq(sk_ranges);
  Index_Partition<Range::POS> sk_pos(sk_ranges);
  // Update strict components
  if (nnc_split)
    points_become_closure_points(sk_eq, sk_strict);
  remove_rows(sk_neg, sk_strict, ns_strict, sat_c_strict);
  // Update nonstrict components
  remove_rows(sk_pos, sk_nonstrict, ns_nonstrict, sat_c_nonstrict);

  return { res_strict, res_nonstrict };
}

} // namespace

/**
 * 2-way rational split method.
 * Letting ph_in be the input value of ph_then, it updates
 *   ph_then = ph_in.add_con(beta)
 *   ph_else = ph_in.add_con(beta_complement)
 * Assumptions:
 *  - ph_then is non-empty and minimized;
 *  - beta is an inequality constraint;
 *  - if ph_then.topology() is CLOSED, then t is CLOSED too
 *    and beta is a non-strict inequality.
 */
void
Q_split_poly(Poly_Impl& ph_then, Poly_Impl& ph_else,
             const Con& beta, Topol t) {
  assert(ph_then.space_dim() == ph_else.space_dim());
  assert(beta.space_dim() <= ph_then.space_dim());
  assert(ph_then.is_minimized());
  assert(ph_else.is_empty());
  const bool nnc_poly = (ph_then.topology() == Topol::NNC);
  assert(t == Topol::CLOSED || nnc_poly);

  if (beta.is_inconsistent()) {
    using std::swap;
    swap(ph_then, ph_else);
    return;
  }
  if (beta.is_tautological() ||
      ph_then.is_empty() ||
      ph_then.space_dim() == 0)
    return;

  ph_else.set_status(Poly_Impl::Status::MINIMIZED);

  const dim_type space_dim = ph_then.space_dim();
  const dim_type num_eq = num_rows(ph_then.cs.sing_rows);

  assert(ph_else.gs.empty());
  // Copy original constraint system
  ph_else.cs = ph_then.cs;

  Con beta_else = complement_con(beta, t);
  const bool strict = beta.is_strict_inequality();
  const bool then_is_strict = (strict || t == Topol::CLOSED);

  // [non]strict suffixes are "real" only if t == Topol::NNC.
  auto& ph_strict = then_is_strict ? ph_then : ph_else;
  auto& ph_nonstrict = then_is_strict ? ph_else : ph_then;
  const auto& beta_strict = then_is_strict ? beta : beta_else;
  const auto& beta_nonstrict = then_is_strict ? beta_else : beta;
  auto& gens_strict = ph_strict.gs;
  auto& gens_nonstrict = ph_nonstrict.gs;
  auto& sat_c_strict = ph_strict.sat_c;
  auto& sat_c_nonstrict = ph_nonstrict.sat_c;
  auto& sat_g_strict = ph_strict.sat_g;
  auto& sat_g_nonstrict = ph_nonstrict.sat_g;
  auto& cons_strict = ph_strict.cs;
  auto& cons_nonstrict = ph_nonstrict.cs;

  // The "strict" side holds the original system.
  if (!then_is_strict) {
    using std::swap;
    swap(gens_strict, gens_nonstrict);
    swap(sat_c_strict, sat_c_nonstrict);
    swap(sat_g_strict, sat_g_nonstrict);
  }

  assert((t == Topol::NNC) == beta_strict.is_strict_inequality());
  const auto [res_strict, res_nonstrict]
    = Q_split_gen_sys(space_dim, num_eq,
                      gens_strict, gens_nonstrict,
                      sat_c_strict, sat_c_nonstrict,
                      sat_g_strict,
                      beta_strict, nnc_poly);
  // For rational splits, these are the only results allowed.
  assert((res_strict == Conv_Iter_Result::OK &&
          res_nonstrict == Conv_Iter_Result::OK)
         ||
         (res_strict == Conv_Iter_Result::EMPTY &&
          res_nonstrict == Conv_Iter_Result::REDUNDANT)
         ||
         (res_strict == Conv_Iter_Result::REDUNDANT &&
          res_nonstrict == Conv_Iter_Result::EMPTY));

  // If the strict side is empty, then `beta_nonstrict' is redundant,
  // so that ph_nonstrict is the same of input polyhedron ph_strict:
  // just swap them.
  if (res_strict == Conv_Iter_Result::EMPTY) {
    using std::swap;
    swap(ph_strict, ph_nonstrict);
    ph_strict.set_empty();
    return;
  }

  // Add the constraints (if not redundant).
  if (res_strict == Conv_Iter_Result::OK) {
    assert(res_nonstrict == Conv_Iter_Result::OK);
    cons_strict.add(beta_strict);
    cons_nonstrict.add(beta_nonstrict);
  }
  // Final simplification.
  if (res_strict == Conv_Iter_Result::EMPTY)
    ph_strict.set_empty();
  else
    simplify(space_dim, cons_strict, sat_c_strict, sat_g_strict,
             gens_strict);
  if (res_nonstrict == Conv_Iter_Result::EMPTY)
    ph_nonstrict.set_empty();
  else
    simplify(space_dim, cons_nonstrict, sat_c_nonstrict, sat_g_nonstrict,
             gens_nonstrict);
}

/******************** Integral split code ********************/

namespace {

void
Z_split_combine(Ranges& sk_ranges_then, Gens& sk_then,
                Integers& sp_then, Sat& sat_then,
                Ranges& sk_ranges_else, Gens& sk_else,
                Integers& sp_else, Sat& sat_else,
                const dim_type min_sat, const dim_type max_sat) {
  assert(sk_ranges_then.size() == sk_ranges_else.size());
  assert(sk_ranges_then.size() == sk_then.size());
  const dim_type old_num_sk_rows = sk_then.size();

  // Note: since initially sk and sat are identical (and they only change
  // by push_back operations), any read to the rows in the range
  // [0, old_num_sk_rows) can be done using sk_then and sat_then.
  assert(sk_then == sk_else);
  assert(sat_then == sat_else);

  // The counting of ones in quick_adj_test() is going to be invoked
  // a quadratic number of times; we linearize and cache it here,
  // before entering the nested loops.
  static PPLITE_TLS Dims sat_ones;
  // Clear and recompute.
  sat_ones.clear();
  for (auto i : index_range(sk_ranges_then)) {
    if (sk_ranges_then[i] == Range::EQ &&
        sk_ranges_else[i] == Range::EQ)
      sat_ones.push_back(0);
    else
      sat_ones.push_back(sat_then[i].count_ones());
  }

  auto get_bits = [](const Ranges& ranges) {
    Bits pos, neg;
    for (auto i : bwd_index_range(ranges)) {
      switch (ranges[i]) {
      case Range::POS:
        pos.set(i);
        break;
      case Range::NEG:
        neg.set(i);
        break;
      case Range::EQ:
      case Range::POS_NEG:
        break;
      }
    }
    return std::make_pair(std::move(pos), std::move(neg));
  };

  auto [pos_then, neg_then] = get_bits(sk_ranges_then);
  auto [pos_else, neg_else] = get_bits(sk_ranges_else);

  // When computing unions, take into account opposite signs in ranges.
  Bits sk_pos = pos_then;
  sk_pos |= neg_else;
  Bits sk_neg = neg_then;
  sk_neg |= pos_else;

  // Note: `pos' and `neg' are relative to the `then' branch.
  for (auto pos : sk_pos) {
    for (auto neg : sk_neg) {
      // Filter out unwanted pos/neg combinations
      if (pos == neg)
        continue;
      const bool combine_then = pos_then[pos] && neg_then[neg];
      // Take into account opposite signs in ranges.
      const bool combine_else = pos_else[neg] && neg_else[pos];
      if (not combine_then && not combine_else)
        continue;

      // Exploit saturation information to perform adjacency tests.
      // Note: delay creation of `new_row' and `new_satrow'.
      const dim_type new_satrow_ones
        = sat_then[pos].count_ones_in_union(sat_then[neg]);
      if (quick_non_adj_test(min_sat, max_sat, new_satrow_ones))
        continue;
      Bits new_satrow = Bits::get_union(sat_then[pos], sat_then[neg]);
      if (not quick_adj_test(pos, neg, new_satrow_ones, sat_ones) &&
          not combinatorial_adj_test(0, old_num_sk_rows,
                                     pos, neg, new_satrow, sat_then))
        continue;

      // Here `pos' and `neg' are adjacent: combine them.
      // New row type is independent of then/else branch.
      auto new_row_type = combine_type(sk_then[pos].type(),
                                       sk_then[neg].type(),
                                       false);
      if (combine_then) {
        if (combine_else)
          // Have to take a copy
          sat_then.add_row(Sat::Row(new_satrow));
        else
          sat_then.add_row(std::move(new_satrow));
        auto new_row = sk_then[neg];
        combine_row(new_row, sk_then[pos], sp_then[neg], sp_then[pos]);
        new_row.set_type(new_row_type);
        sk_then.push_back(std::move(new_row));
      }
      if (combine_else) {
        sat_else.add_row(std::move(new_satrow));
        // Take into account opposite ranges.
        auto new_row = sk_else[pos];
        combine_row(new_row, sk_else[neg], sp_else[pos], sp_else[neg]);
        new_row.set_type(new_row_type);
        sk_else.push_back(std::move(new_row));
      }
    } // end loop on sk_neg
  } // end loop on sk_pos

  sk_ranges_then.resize(sk_then.size(), Range::EQ);
  sk_ranges_else.resize(sk_else.size(), Range::EQ);
}

// Incrementally adjust scalar products based on delta_inhomo.
inline void
Z_adjust_scal_prods(Integers& sps, const Gens& sk,
                    const Integer& delta_inhomo) {
  for (auto i : index_range(sps)) {
    neg_assign(sps[i]);
    if (sk[i].is_point())
      add_mul_assign(sps[i], delta_inhomo, sk[i].divisor());
  }
}

inline void
Z_adjust_sat(Ranges& sk_ranges, Sat& sat) {
  sat.add_cols(1);
  Index_Partition<Range::POS> sk_pos(sk_ranges);
  update_last_sat_column(sat, sk_pos);
}

inline void
Z_cleanup(Ranges& sk_ranges, Gens& sk, Sat& sat) {
  NS_Rows ns_dummy;
  Index_Partition<Range::NEG> sk_neg(sk_ranges);
  remove_rows(sk_neg, sk, ns_dummy, sat);
}

inline void
Z_adjust_and_cleanup(const Con& beta,
                     Ranges& sk_ranges, Gens& sk, Sat& sat) {
  if (beta.is_equality()) {
    // equality: sat should not be adjusted;
    // modify sk_ranges so that rows in SK+ will get removed.
    std::replace(sk_ranges.begin(), sk_ranges.end(),
                 Range::POS, Range::NEG);
  } else {
    // inequality: sat should be adjusted.
    Z_adjust_sat(sk_ranges, sat);
  }
  Z_cleanup(sk_ranges, sk, sat);
}

inline void
Z_single_combine(const Con& beta, dim_type min_sat, dim_type max_sat,
                 Ranges& sk_ranges, Gens& sk, Integers& sp, Sat& sat) {
  combine_sk_rows(sk_ranges, sk, sp, sat, false, min_sat, max_sat);
  Z_adjust_and_cleanup(beta, sk_ranges, sk, sat);
}

/**
 * Computes the split sk systems, copying the original one.
 * NB: sk/sat_[then|else] are consistent only if result is not EMPTY.
 * If inconsistent, it is up to the caller to set the empty DD.
 */
Split_Result
Z_split_skel_on_ineq(dim_type space_dim, dim_type num_lines, dim_type num_eq,
                     Gens& sk_then, Gens& sk_else,
                     Sat& sat_then, Sat& sat_else,
                     const Con& beta_then, const Con& beta_else) {
  assert(beta_then.is_nonstrict_inequality() &&
         beta_else.is_nonstrict_inequality());
  // Note: sk_then and sat_then are in-out params;
  // sk_else and sat_else are out params.
  // We want to be lazy in populating the *_else params,
  // so as to factor common computations as much as possible.
  assert(sk_else.empty() && sat_else.empty());

  // Start processing the `then' component.
  Integers scal_prods_then;
  sp::compute_prods(beta_then, sk_then, scal_prods_then);

  Ranges sk_ranges_then;
  sk_partition(scal_prods_then, sk_ranges_then);
  Index_Partition<Range::NEG> sk_neg_then(sk_ranges_then);
  Index_Partition<Range::EQ> sk_eq_then(sk_ranges_then);
  Index_Partition<Range::POS> sk_pos_then(sk_ranges_then);
  bool combine_then = not sk_neg_then.empty() && not sk_pos_then.empty();
  // Check for emptiness/redundancy.
  auto res_then = combine_then ? Conv_Iter_Result::OK
    : sk_neg_then.empty() ? Conv_Iter_Result::REDUNDANT
    : (sk_pos_then.empty() && has_no_eq_point(sk_ranges_then, sk_then)
       ? Conv_Iter_Result::EMPTY
       : Conv_Iter_Result::OK);

  Integers scal_prods_else;
  if (combine_then) {
    // Need deep copy.
    scal_prods_else = scal_prods_then;
  } else {
    // Move (no need to deep copy)
    scal_prods_else = std::move(scal_prods_then);
  }
  Integer delta_inhomo = beta_then.inhomo_term() + beta_else.inhomo_term();
  Z_adjust_scal_prods(scal_prods_else, sk_then, delta_inhomo);

  // Now process `else' component, avoiding duplication of work.
  Ranges sk_ranges_else;
  sk_partition(scal_prods_else, sk_ranges_else);
  Index_Partition<Range::NEG> sk_neg_else(sk_ranges_else);
  Index_Partition<Range::EQ> sk_eq_else(sk_ranges_else);
  Index_Partition<Range::POS> sk_pos_else(sk_ranges_else);
  bool combine_else = not sk_neg_else.empty() && not sk_pos_else.empty();

  auto res_else = combine_else ? Conv_Iter_Result::OK
    : sk_neg_else.empty() ? Conv_Iter_Result::REDUNDANT
    : (sk_pos_else.empty() && has_no_eq_point(sk_ranges_else, sk_then)
       ? Conv_Iter_Result::EMPTY
       : Conv_Iter_Result::OK);

  // Eager return in case both are empty.
  if (res_then == Conv_Iter_Result::EMPTY &&
      res_else == Conv_Iter_Result::EMPTY)
    return { res_then, res_else };

  // Needed when combining.
  const dim_type min_sat = space_dim - num_lines - 1;
  const dim_type max_sat = sat_then.num_cols() + num_eq;

  if (combine_then && combine_else) {
    // Combine both.
    sk_else = sk_then;
    sat_else = sat_then;
    Z_split_combine(sk_ranges_then, sk_then, scal_prods_then, sat_then,
                    sk_ranges_else, sk_else, scal_prods_else, sat_else,
                    min_sat, max_sat);
    Z_adjust_and_cleanup(beta_then, sk_ranges_then, sk_then, sat_then);
    Z_adjust_and_cleanup(beta_else, sk_ranges_else, sk_else, sat_else);
    return { res_then, res_else };
  }

  if (combine_then) {
    assert(not combine_else);
    // To be computed before combining `then'
    if (res_else != Conv_Iter_Result::EMPTY) {
      sk_else = sk_then;
      sat_else = sat_then;
      if (res_else == Conv_Iter_Result::OK)
        Z_adjust_and_cleanup(beta_else, sk_ranges_else, sk_else, sat_else);
    }
    // Now combine `then'
    Z_single_combine(beta_then, min_sat, max_sat,
                     sk_ranges_then, sk_then, scal_prods_then, sat_then);
    return { res_then, res_else };
  }

  if (combine_else) {
    assert(not combine_then);
    // To be computed before combining `else'
    if (res_then == Conv_Iter_Result::EMPTY) {
      // can avoid copy
      sk_else = std::move(sk_then);
      sat_else = std::move(sat_then);
    } else {
      // perform copy
      sk_else = sk_then;
      sat_else = sat_then;
      if (res_then == Conv_Iter_Result::OK)
        Z_adjust_and_cleanup(beta_then, sk_ranges_then, sk_then, sat_then);
    }
    // Now combine `else'
    Z_single_combine(beta_else, min_sat, max_sat,
                     sk_ranges_else, sk_else, scal_prods_else, sat_else);
    return { res_then, res_else };
  }

  // No combine at all: they can't be both empty.
  if (res_else != Conv_Iter_Result::EMPTY) {
    if (res_then == Conv_Iter_Result::EMPTY) {
      // avoid deep copy
      sk_else = std::move(sk_then);
      sat_else = std::move(sat_then);
    } else {
      // need to copy
      sk_else = sk_then;
      sat_else = sat_then;
    }
  }
  if (res_then == Conv_Iter_Result::OK)
    Z_adjust_and_cleanup(beta_then, sk_ranges_then, sk_then, sat_then);
  if (res_else == Conv_Iter_Result::OK)
    Z_adjust_and_cleanup(beta_else, sk_ranges_else, sk_else, sat_else);
  return { res_then, res_else };
}

Split_Result
Z_split_skel_on_eq(dim_type space_dim, dim_type num_lines, dim_type num_eq,
                   Gens& sk_equal, Gens& sk_ineq,
                   Sat& sat_equal, Sat& sat_ineq,
                   const Con& beta_equal, Con& beta_ineq) {
  assert(beta_equal.is_equality());
  // Note: sk_equal and sat_equal are in-out params;
  // sk_ineq, sat_ineq and beta_ineq are out params.
  // We want to be lazy in populating the *_ineq params,
  // so as to factor common computations as much as possible.
  assert(sk_ineq.empty() && sat_ineq.empty());

  // Start processing the `equal' (i.e., the `then' branch) component.
  Integers scal_prods_equal;
  sp::compute_prods(beta_equal, sk_equal, scal_prods_equal);

  // Check for special cases.
  const int sp_sign = sgn(scal_prods_equal[0]);
  const bool same_sign
    = all_of(scal_prods_equal,
             [sp_sign](const Integer& sp) { return sgn(sp) == sp_sign; });
  if (same_sign) {
    // All products have the same sign:
    // the equality is either redundant or inconsistent
    if (sp_sign == 0) {
      // res_equal is included in beta_equal
      return { Conv_Iter_Result::REDUNDANT, Conv_Iter_Result::EMPTY };
    } else {
      // res_equal is disjoint from beta_equal:
      // move equal poly into else poly
      sk_ineq = std::move(sk_equal);
      sat_ineq = std::move(sat_equal);
      return { Conv_Iter_Result::EMPTY, Conv_Iter_Result::REDUNDANT };
    }
  }

  Ranges sk_ranges_equal;
  sk_partition(scal_prods_equal, sk_ranges_equal);
  Index_Partition<Range::NEG> sk_neg_equal(sk_ranges_equal);
  Index_Partition<Range::EQ> sk_eq_equal(sk_ranges_equal);
  Index_Partition<Range::POS> sk_pos_equal(sk_ranges_equal);
  bool combine_equal = not sk_neg_equal.empty() && not sk_pos_equal.empty();
  // Check for emptiness/redundancy.
  auto res_equal = combine_equal ? Conv_Iter_Result::OK
    : sk_neg_equal.empty() ? Conv_Iter_Result::REDUNDANT
    : (sk_pos_equal.empty() && has_no_eq_point(sk_ranges_equal, sk_equal)
       ? Conv_Iter_Result::EMPTY
       : Conv_Iter_Result::OK);

  // Needed when combining.
  const dim_type min_sat = space_dim - num_lines - 1;
  const dim_type max_sat = sat_equal.num_cols() + num_eq;

  Integers scal_prods_ineq;

  // Check for lucky cases.
  auto [beta_lt, beta_gt] = integral_complement_eq(beta_equal);
  // Check if ph_lt happens to be empty (lucky case).
  Integers scal_prods_lt = scal_prods_equal;
  Integer delta_inhomo = beta_equal.inhomo_term() + beta_lt.inhomo_term();
  Z_adjust_scal_prods(scal_prods_lt, sk_equal, delta_inhomo);
  const bool ph_lt_empty
    = all_of(scal_prods_lt, [](const Integer& sp) { return sgn(sp) < 0; });
  if (ph_lt_empty) {
    // ph_lt is empty: the else branch can be based on beta_gt.
    beta_ineq = beta_gt;
    scal_prods_ineq = scal_prods_lt;
    delta_inhomo = beta_lt.inhomo_term() + beta_gt.inhomo_term();
    Z_adjust_scal_prods(scal_prods_ineq, sk_equal, delta_inhomo);
  } else {
    // Check if ph_gt happens to be empty (other lucky case).
    Integers scal_prods_gt = scal_prods_lt;
    delta_inhomo = beta_lt.inhomo_term() + beta_gt.inhomo_term();
    Z_adjust_scal_prods(scal_prods_gt, sk_equal, delta_inhomo);
    const bool ph_gt_empty
      = all_of(scal_prods_gt, [](const Integer& sp) { return sgn(sp) < 0; });
    if (ph_gt_empty) {
      // ph_gt is empty: the else branch can be based on beta_lt.
      beta_ineq = beta_lt;
      scal_prods_ineq = std::move(scal_prods_lt);
    } else {
      // unlucky case: ph_lt and ph_gt are both not empty,
      // so that we cannot filter on the else branch.
      auto res_ineq = Conv_Iter_Result::REDUNDANT;
      sk_ineq = sk_equal;
      sat_ineq = sat_equal;
      if (combine_equal) {
        Z_single_combine(beta_equal, min_sat, max_sat, sk_ranges_equal,
                         sk_equal, scal_prods_equal, sat_equal);
      }
      return { res_equal, res_ineq };
    }
  }

  // Now process `else' component, avoiding duplication of work.
  Ranges sk_ranges_ineq;
  sk_partition(scal_prods_ineq, sk_ranges_ineq);
  Index_Partition<Range::NEG> sk_neg_ineq(sk_ranges_ineq);
  Index_Partition<Range::EQ> sk_eq_ineq(sk_ranges_ineq);
  Index_Partition<Range::POS> sk_pos_ineq(sk_ranges_ineq);
  bool combine_ineq = not sk_neg_ineq.empty() && not sk_pos_ineq.empty();

  auto res_ineq = combine_ineq ? Conv_Iter_Result::OK
    : sk_neg_ineq.empty() ? Conv_Iter_Result::REDUNDANT
    : (sk_pos_ineq.empty() && has_no_eq_point(sk_ranges_ineq, sk_equal)
       ? Conv_Iter_Result::EMPTY
       : Conv_Iter_Result::OK);

  // Eager return in case both are empty.
  if (res_equal == Conv_Iter_Result::EMPTY &&
      res_ineq == Conv_Iter_Result::EMPTY)
    return { res_equal, res_ineq };

  if (combine_equal && combine_ineq) {
    // Deep copy and then combine both.
    sk_ineq = sk_equal;
    sat_ineq = sat_equal;
    Z_split_combine(sk_ranges_equal, sk_equal, scal_prods_equal, sat_equal,
                    sk_ranges_ineq, sk_ineq, scal_prods_ineq, sat_ineq,
                    min_sat, max_sat);
    Z_adjust_and_cleanup(beta_equal, sk_ranges_equal, sk_equal, sat_equal);
    Z_adjust_and_cleanup(beta_ineq, sk_ranges_ineq, sk_ineq, sat_ineq);
    return { res_equal, res_ineq };
  }

  if (combine_equal) {
    assert(not combine_ineq);
    // To be computed before combining `equal'
    if (res_ineq != Conv_Iter_Result::EMPTY) {
      sk_ineq = sk_equal;
      sat_ineq = sat_equal;
      if (res_ineq == Conv_Iter_Result::OK)
        Z_adjust_and_cleanup(beta_ineq, sk_ranges_ineq, sk_ineq, sat_ineq);
    }
    // Now combine `equal'.
    Z_single_combine(beta_equal, min_sat, max_sat, sk_ranges_equal,
                     sk_equal, scal_prods_equal, sat_equal);
    return { res_equal, res_ineq };
  }

  if (combine_ineq) {
    assert(not combine_equal);
    // To be computed before combining `else'
    if (res_equal == Conv_Iter_Result::EMPTY) {
      // avoid deep copy
      sk_ineq = std::move(sk_equal);
      sat_ineq = std::move(sat_equal);
    } else {
      // deep copy
      sk_ineq = sk_equal;
      sat_ineq = sat_equal;
      if (res_equal == Conv_Iter_Result::OK)
        Z_adjust_and_cleanup(beta_equal, sk_ranges_equal, sk_equal, sat_equal);
    }
    // Now combine `ineq'
    Z_single_combine(beta_ineq, min_sat, max_sat,
                     sk_ranges_ineq, sk_ineq, scal_prods_ineq, sat_ineq);
    return { res_equal, res_ineq };
  }

  // No combine at all: they can't be both empty.
  if (res_ineq != Conv_Iter_Result::EMPTY) {
    if (res_equal == Conv_Iter_Result::EMPTY) {
      // avoid deep copy
      sk_ineq = std::move(sk_equal);
      sat_ineq = std::move(sat_equal);
    } else {
      // deep copy
      sk_ineq = sk_equal;
      sat_ineq = sat_equal;
    }
  }
  if (res_equal == Conv_Iter_Result::OK)
    Z_adjust_and_cleanup(beta_equal, sk_ranges_equal, sk_equal, sat_equal);
  if (res_ineq == Conv_Iter_Result::OK)
    Z_adjust_and_cleanup(beta_ineq, sk_ranges_ineq, sk_ineq, sat_ineq);
  return { res_equal, res_ineq };
}

Split_Result
Z_split_gen_sys(dim_type space_dim, dim_type num_eq,
                Gen_Sys& gen_sys_then, Gen_Sys& gen_sys_else,
                Sat& sat_c_then, Sat& sat_c_else,
                const Con& beta_then, Con& beta_else) {
  Gens& sg_then = gen_sys_then.sing_rows;
  Gens& sg_else = gen_sys_else.sing_rows;
  const bool split_on_eq = beta_then.is_equality();

  // Check if all singular rows saturate beta_then.
  const dim_type ex = check_sing_rows(beta_then, sg_then);
  if (ex < num_rows(sg_then)) {
    // NOTE: optimizable; here we just copy and work twice.
    gen_sys_else = gen_sys_then;
    sat_c_else = sat_c_then;
    Ranges dummy;
    process_violating_singular(ex, beta_then,
                               gen_sys_then, sat_c_then,
                               dummy, dummy);
    if (split_on_eq) {
      // Here beta_else is a dummy constraint.
      // Since a line is violating equality beta_then,
      // the else branch cannot be refined;
      // hence we flag the (dummy) beta_else constraint as redundant.
      return { Conv_Iter_Result::OK, Conv_Iter_Result::REDUNDANT };
    } else {
      // Here beta_else is not a dummy constraint.
      process_violating_singular(ex, beta_else,
                                 gen_sys_else, sat_c_else,
                                 dummy, dummy);
      return { Conv_Iter_Result::OK, Conv_Iter_Result::OK };
    }
  }

  // Here all lines saturate the constraint.
  // Copy lines
  sg_else = sg_then;
  dim_type num_lines = num_rows(sg_then);

  // Split the skeleton part.
  Gens& sk_then = gen_sys_then.sk_rows;
  Gens& sk_else = gen_sys_else.sk_rows;
  return split_on_eq
    ? Z_split_skel_on_eq(space_dim, num_lines, num_eq,
                         sk_then, sk_else,
                         sat_c_then, sat_c_else,
                         beta_then, beta_else)
    : Z_split_skel_on_ineq(space_dim, num_lines, num_eq,
                           sk_then, sk_else,
                           sat_c_then, sat_c_else,
                           beta_then, beta_else);
}

void
Z_split_simplify(Conv_Iter_Result res, Poly_Impl& ph, const Con& beta) {
  switch (res) {
  case Conv_Iter_Result::OK:
    ph.cs.add(beta);
    simplify(ph.space_dim(), ph.cs, ph.sat_c, ph.sat_g, ph.gs);
    break;
  case Conv_Iter_Result::EMPTY:
    ph.set_empty();
    break;
  case Conv_Iter_Result::REDUNDANT:
    // nothing to do
    break;
  }
}

void
Z_split_simplify_both(Split_Result res,
                      Poly_Impl& ph_then, const Con& beta_then,
                      Poly_Impl& ph_else, const Con& beta_else) {
  const auto [res_then, res_else] = res;

  // If needed, copy/move cons_then into cons_else.
  if (res_else != Conv_Iter_Result::EMPTY) {
    if (res_then == Conv_Iter_Result::EMPTY) {
      ph_else.cs = std::move(ph_then.cs);
      ph_else.sat_g = std::move(ph_then.sat_g);
    } else {
      ph_else.cs = ph_then.cs;
      ph_else.sat_g = ph_else.sat_c.transpose();
    }
  }
  Z_split_simplify(res_then, ph_then, beta_then);
  Z_split_simplify(res_else, ph_else, beta_else);
}

/**
 * 2-way integral split method on an inequality constraint.
 * Letting
 *  - ph_in be the input value of ph_then,
 *  - beta_then be the integral refinement of beta,
 *  - beta_else be the integral complement of beta_then
 * it updates
 *   ph_then = ph_in.add_con(beta_then)
 *   ph_else = ph_in.add_con(beta_else)
 * Assumptions:
 *  - ph_then.topology() is CLOSED
 *  - ph_then is non-empty and minimized;
 *  - beta is a (strict or non-strict) inequality constraint;
 */
void
Z_split_poly_on_ineq(Poly_Impl& ph_then, Poly_Impl& ph_else,
                     const Con& beta) {
  assert(ph_then.is_necessarily_closed());
  assert(ph_else.is_necessarily_closed());
  assert(ph_then.space_dim() == ph_else.space_dim());
  assert(beta.space_dim() <= ph_then.space_dim());
  assert(not beta.is_equality());
  assert(ph_then.is_minimized() && not ph_then.is_empty());
  assert(ph_else.is_empty());

  if (beta.is_inconsistent()) {
    using std::swap;
    swap(ph_then, ph_else);
    return;
  }
  if (beta.is_tautological() ||
      ph_then.is_empty() ||
      ph_then.space_dim() == 0)
    return;

  auto [beta_then, beta_else] = integral_complement_cons(beta);
  assert(beta_then.is_nonstrict_inequality());
  assert(beta_else.is_nonstrict_inequality());

  const dim_type space_dim = ph_then.space_dim();
  const dim_type num_eq = num_rows(ph_then.cs.sing_rows);
  ph_else.set_status(Poly_Impl::Status::MINIMIZED);
  // Note: in the integral case both results may be empty.
  auto split_res = Z_split_gen_sys(space_dim, num_eq,
                                   ph_then.gs, ph_else.gs,
                                   ph_then.sat_c, ph_else.sat_c,
                                   beta_then, beta_else);
  Z_split_simplify_both(split_res,
                        ph_then, beta_then,
                        ph_else, beta_else);
}

/**
 * 2-way integral split on an equality constraint.
 * Letting
 *  - ph_in be the input value of ph_eq,
 *  - beta_eq := (expr == 0) be an equality constraint,
 *  - beta_lt be the integral refinement of (expr < 0),
 *  - beta_gt be the integral refinement of (expr > 0)
 * it updates
 *   ph_eq = ph_in.add_con(beta_eq)
 *   ph_ineq = ph_in.add_con(beta_lt), if ph_in and beta_gt are disjoint;
 *           = ph_in.add_con(beta_gt), if ph_in and beta_lt are disjoint;
 *           = ph_in,                  otherwise.
 * Assumptions:
 *  - ph_eq.topology() is CLOSED
 *  - ph_eq is non-empty and minimized;
 *  - ph_ineq is empty, CLOSED and same space dim as ph_eq
 *  - beta_eq is an equality constraint.
 */
void
Z_split_poly_on_eq(Poly_Impl& ph_eq, Poly_Impl& ph_ineq,
                   const Con& beta_eq) {
  assert(ph_eq.is_necessarily_closed());
  assert(ph_ineq.is_necessarily_closed());
  assert(ph_eq.space_dim() == ph_ineq.space_dim());
  assert(ph_eq.is_minimized() && not ph_eq.is_empty());
  assert(ph_ineq.is_empty());
  assert(beta_eq.is_equality());
  assert(beta_eq.space_dim() <= ph_eq.space_dim());

  if (is_integral_inconsistent(beta_eq)) {
    using std::swap;
    swap(ph_eq, ph_ineq);
    return;
  }
  if (beta_eq.is_tautological() ||
      ph_eq.is_empty() ||
      ph_eq.space_dim() == 0)
    return;

  const dim_type space_dim = ph_eq.space_dim();
  const dim_type num_eq = num_rows(ph_eq.cs.sing_rows);
  ph_ineq.set_status(Poly_Impl::Status::MINIMIZED);
  // dummy input: it will be overwritten
  auto beta_ineq = Con::zero_dim_false();
  // Note: in the integral case with an equality,
  // the result OK-REDUNDANT is allowed.
  auto split_res = Z_split_gen_sys(space_dim, num_eq,
                                   ph_eq.gs, ph_ineq.gs,
                                   ph_eq.sat_c, ph_ineq.sat_c,
                                   beta_eq, beta_ineq);
  Z_split_simplify_both(split_res, ph_eq, beta_eq, ph_ineq, beta_ineq);
}

} // namespace

void
Z_split_poly(Poly_Impl& ph_then, Poly_Impl& ph_else,
             const Con& beta) {
  if (beta.is_equality())
    Z_split_poly_on_eq(ph_then, ph_else, beta);
  else
    Z_split_poly_on_ineq(ph_then, ph_else, beta);
}

} // namespace detail
} // namespace pplite
