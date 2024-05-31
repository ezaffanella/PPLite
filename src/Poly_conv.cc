/* PPLite: a lightweight library for convex polyhedra derived from PPL.
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

#define PPLITE_SORT_PENDING 0

#define PPLITE_FAST_PATH_STATS 0

#include "globals.hh"
#include "support_utils.hh"
#include "utils.hh"
#include "Bits.hh"
#include "Con.hh"
#include "Gen.hh"
#include "Index_Partition.hh"
#include "Integer.hh"
#include "Poly_min.hh"
#include "Sat.hh"
#include "Scalar_Prod.hh"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <numeric> // for std::iota
#include <utility>
#include <vector>

namespace pplite {
namespace detail {

#if PPLITE_FAST_PATH_STATS

#define FAST_PATH_STATS_INCR_CALLS (++fp_tests)
#define FAST_PATH_STATS_INCR_TRUE (++fp_tests_true)
#define FAST_PATH_STATS_INCR_DOUBT (++fp_tests_doubt)

#define FAST_PATH_STATS_INCR_EXEC_FAST \
  do { if (!con_to_gen) ++fp_exec_fast; } while (false)
#define FAST_PATH_STATS_INCR_EXEC_SLOW \
  do { if (!con_to_gen) ++fp_exec_slow; } while (false)

#define FAST_PATH_STATS_RESET \
  do { if (!con_to_gen) fp_stats_reset(); } while (false)
#define FAST_PATH_STATS_PRINT \
  do { if (!con_to_gen) fp_stats_print(); } while (false)

namespace {

PPLITE_TLS long fp_tests = 0;
PPLITE_TLS long fp_tests_true = 0;
PPLITE_TLS long fp_tests_doubt = 0;
PPLITE_TLS long fp_exec_fast = 0;
PPLITE_TLS long fp_exec_slow = 0;

inline void fp_stats_reset() {
  fp_tests = 0;
  fp_tests_true = 0;
  fp_tests_doubt = 0;
  fp_exec_fast = 0;
  fp_exec_slow = 0;
}

inline void fp_stats_print() {
  auto get_perc = [](unsigned num) -> long {
    return lround((num / static_cast<double>(fp_tests)) * 100.0);
  };
  if (fp_tests == 0)
    return;
  auto& os = std::cerr;
  auto perc_true = get_perc(fp_tests_true);
  auto perc_doubt = get_perc(fp_tests_doubt);
  auto perc_fast =  get_perc(fp_exec_fast);
  auto perc_slow = get_perc(fp_exec_slow);
  os << "\n"
     << "FAST PATH:"
     << " calls " << fp_tests
     << " --- "
     << " %true " << perc_true
     << " %doubt " << perc_doubt
     << " --- "
     << " %fast " << perc_fast
     << " %slow " << perc_slow
     << "\n";
}

} // namespace

#else // !PPLITE_FAST_PATH_STATS

#define FAST_PATH_STATS_INCR_CALLS
#define FAST_PATH_STATS_INCR_TRUE
#define FAST_PATH_STATS_INCR_DOUBT
#define FAST_PATH_STATS_INCR_EXEC_FAST
#define FAST_PATH_STATS_INCR_EXEC_SLOW
#define FAST_PATH_STATS_RESET
#define FAST_PATH_STATS_PRINT

#endif // !PPLITE_FAST_PATH_STATS


#if PPLITE_NOISY_CONVERSION

PPLITE_TLS unsigned long num_quick_nonadj = 0;
PPLITE_TLS unsigned long num_quick_adj = 0;
PPLITE_TLS unsigned long num_comb_nonadj = 0;
PPLITE_TLS unsigned long num_comb_adj = 0;

#define INCR_QUICK_NONADJ (++num_quick_nonadj)
#define INCR_QUICK_ADJ    (++num_quick_adj)
#define INCR_COMB_NONADJ  (++num_comb_nonadj)
#define INCR_COMB_ADJ     (++num_comb_adj)

template <typename Dst_Sys>
void
output_counters(const char* kind, dim_type k, const Dst_Sys& dst) {
  auto sg_count = dst.sing_rows.size();
  auto sk_count = dst.sk_rows.size();
  auto ns_count = dst.ns_rows.size();

  std::ostream& os = std::cerr;
  const auto num_adj
    = num_quick_nonadj + num_quick_adj + num_comb_nonadj + num_comb_adj;
  os << "  adj tests = " << num_adj;
  os << " (quick no = " << num_quick_nonadj
     << ", quick yes = " << num_quick_adj
     << ", comb no = " << num_comb_nonadj
     << ", comb yes = " << num_comb_adj
     << ")\n\n";
  // Reset counters.
  num_quick_nonadj = 0;
  num_quick_adj = 0;
  num_comb_nonadj = 0;
  num_comb_adj = 0;
  os << "conv iter " << kind << " k = " << k
     << " | sg = " << sg_count
     << " | sk = " << sk_count
     << " | ns = " << ns_count
     << "\n";
}

// Explicit template instantiations.
template void
output_counters<Con_Sys>(const char*, dim_type, const Con_Sys&);
template void
output_counters<Gen_Sys>(const char*, dim_type, const Gen_Sys&);

#else // !PPLITE_NOISY_CONVERSION

#define INCR_QUICK_NONADJ
#define INCR_QUICK_ADJ
#define INCR_COMB_NONADJ
#define INCR_COMB_ADJ

#endif // !PPLITE_NOISY_CONVERSION

// Partitions skeleton rows into ranges (stored in `sk_ranges')
// based on sign of scalar products.
void
sk_partition(const Integers& sp, Ranges& sk_ranges) {
  sk_ranges.resize(sp.size());
  for (auto i : bwd_index_range(sp)) {
    const int sp_sign = sgn(sp[i]);
    sk_ranges[i]
      = (sp_sign == 0) ? Range::EQ
      : (sp_sign == 1) ? Range::POS
      : Range::NEG;
  }
}

void
sk_ranges_to_index_sets(const Ranges& sk_ranges,
                        Index_Set& iset_neg,
                        Index_Set& iset_eq,
                        Index_Set& iset_pos) {
  assert(iset_neg.empty());
  assert(iset_eq.empty());
  assert(iset_pos.empty());
  for (auto i : bwd_index_range(sk_ranges)) {
    switch (sk_ranges[i]) {
    case Range::NEG:
      iset_neg.set(i);
      break;
    case Range::EQ:
      iset_eq.set(i);
      break;
    case Range::POS:
      iset_pos.set(i);
      break;
    default:
      PPLITE_UNREACH;
    }
  }
}

void
sk_ranges_to_noneq_index_sets(const Ranges& sk_ranges,
                              Index_Set& iset_neg,
                              Index_Set& iset_pos) {
  assert(iset_neg.empty());
  assert(iset_pos.empty());
  for (auto i : bwd_index_range(sk_ranges)) {
    switch (sk_ranges[i]) {
    case Range::NEG:
      iset_neg.set(i);
      break;
    case Range::POS:
      iset_pos.set(i);
      break;
    case Range::EQ:
    default:
      break;
    }
  }
}

Index_Set
sk_ranges_to_eq_index_set(const Ranges& sk_ranges) {
  Index_Set res;
  for (auto i : bwd_index_range(sk_ranges))
    if (sk_ranges[i] == Range::EQ)
      res.set(i);
  return res;
}

// Search range `[first, last)' for a generator,
// different from `pos' and `neg',
// saturating all the constraints saturated by both `pos' and `neg':
// if found, `pos' and `neg' are not adjacent.
// Note: the test is symmetric; names `pos' and `neg' are arbitrary.
bool
combinatorial_adj_test(const dim_type first,
                       const dim_type last,
                       const dim_type pos,
                       const dim_type neg,
                       const Bits& new_satrow,
                       const Sat& sat_src) {
  assert(first < last && pos != neg);
  assert(first <= pos && pos < last);
  assert(first <= neg && neg < last);
  const dim_type lower = std::min(pos, neg);
  const dim_type upper = std::max(pos, neg);
  for (auto i : range(first, lower)) {
    if (subset_eq(sat_src[i], new_satrow))
      return false;
  }
  for (auto i : range(lower+1, upper)) {
    if (subset_eq(sat_src[i], new_satrow))
      return false;
  }
  for (auto i : range(upper+1, last)) {
    if (subset_eq(sat_src[i], new_satrow))
      return false;
  }
  return true;
}

// Note: here we do NOT map points into closure points,
// even when the processed constraint is a strict inequality.
// This mapping will be done by the caller.
// Note: this function works for both
//   Src_Row = Con and SK_Rows = Gens
// or
//   Src_Row = Gen and SK_Rows = Cons.
template <typename Src_Row, typename SK_Rows>
void
combine_with_ex_sing(const dim_type ex,
                     const Integer& sp_ex,
                     const Src_Row& src_k,
                     SK_Rows& sg_rows,
                     SK_Rows& sk_rows) {
  assert(0 <= ex && ex < num_rows(sg_rows));
  assert(is_ex_singular(sg_rows[ex]) && sp_ex > 0);
  // Allocating temporaries once for all before entering the loop.
  Integer sp_i;
  const auto& sg_ex = sg_rows[ex];
  Index_Set src_k_nz;
  const bool sparse = sp::is_sparse(src_k, src_k_nz);
  for (auto i : range(ex+1, num_rows(sg_rows))) {
    auto& row_i = sg_rows[i];
    if (sparse)
      sp::assign(sp_i, src_k_nz, src_k, row_i);
    else
      sp::assign(sp_i, src_k, row_i);
    if (sp_i != 0)
      combine_row(row_i, sg_ex, sp_i, sp_ex);
  }
  for (auto& row_i : sk_rows) {
    if (sparse)
      sp::assign(sp_i, src_k_nz, src_k, row_i);
    else
      sp::assign(sp_i, src_k, row_i);
    if (sp_i != 0)
      combine_row(row_i, sg_ex, sp_i, sp_ex);
  }
}

// Combines the skeleton generator in SK+ with those in SK-.
// Note: this works as is also for constraints,
// provided is_strict_ineq actually stores is_point.
template <typename SK_Rows>
void
combine_sk_rows(Ranges& sk_ranges, SK_Rows& sk_rows,
                Integers& sp, Sat& sat_src, const bool is_strict_ineq,
                const dim_type min_sat, const dim_type max_sat) {
  const dim_type old_num_sk_rows = sk_rows.size();

  // Using dense Index_Partition for sk_pos and sk_neg causes slow downs.
  static PPLITE_TLS Dims sk_pos, sk_neg;
  // The counting of ones in quick_adj_test() is going to be invoked
  // a quadratic number of times; we linearize and cache it here,
  // before entering the nested loops.
  static PPLITE_TLS Dims sat_ones;

  // Clear and recompute.
  sk_pos.clear();
  sk_neg.clear();
  sat_ones.clear();
  for (auto i : index_range(sk_ranges)) {
    if (sk_ranges[i] == Range::EQ) {
      sat_ones.push_back(0);
      continue;
    }
    sat_ones.push_back(sat_src[i].count_ones());
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
        = sat_src[pos].count_ones_in_union(sat_src[neg]);
      if (quick_non_adj_test(min_sat, max_sat, new_satrow_ones)) {
        INCR_QUICK_NONADJ;
        continue;
      }
      Bits new_satrow = Bits::get_union(sat_src[pos], sat_src[neg]);
      if (quick_adj_test(pos, neg, new_satrow_ones, sat_ones)) {
        INCR_QUICK_ADJ;
        // Adjacent: fall through.
      } else if (combinatorial_adj_test(0, old_num_sk_rows,
                                        pos, neg, new_satrow, sat_src)) {
        INCR_COMB_ADJ;
        // Adjacent: fall through.
      } else {
        INCR_COMB_NONADJ;
        continue;
      }

      // Here `pos' and `neg' are adjacent: update `sat_src'.
      sat_src.add_row(std::move(new_satrow));
      // Add the new row to `sk_rows'
      auto new_row = sk_rows[neg];
      combine_row(new_row, sk_rows[pos], sp[neg], sp[pos]);
      new_row.set_type(combine_type(sk_rows[pos].type(),
                                    sk_rows[neg].type(),
                                    is_strict_ineq));
      sk_rows.push_back(std::move(new_row));
      assert(num_rows(sk_rows) == sat_src.num_rows());

    } // end loop on sk_neg
  } // end loop on sk_pos

  sp.resize(sk_rows.size());
  sk_ranges.resize(sk_rows.size(), Range::EQ);
}

/* Returns true if new_ns is NOT redundant in set_of_ns */
bool
remove_redundants(const Index_Set& new_ns, NS_Rows& set_of_ns) {
  Dims to_be_del;
  auto i = 0, i_end = num_rows(set_of_ns);
  for ( ; i != i_end; ++i) {
    const auto& ns = set_of_ns[i];
    if (subset_eq(ns, new_ns)) {
      // ns makes new_ns redundant: nothing to do
      assert(to_be_del.empty());
      return false;
    }
    if (subset_eq(new_ns, ns)) {
      assert(new_ns != ns);
      // new_ns makes ns redundant: remove it
      to_be_del.push_back(i);
      break;
    }
  }
  // Here new_ns is not redundant
  if (i != i_end) {
    assert(not to_be_del.empty());
    // Check for other redundancies
    for (++i; i != i_end; ++i) {
      const auto& ns = set_of_ns[i];
      if (subset_eq(new_ns, ns)) {
        assert(new_ns != ns);
        to_be_del.push_back(i);
      }
    }
    erase_using_sorted_indices(set_of_ns, to_be_del);
  }
  return true;
}

// If not redundant, adds `new_ns' in `set_of_ns' (minimizing it).
// We provide the copy/move variants.
void
add_to_minimal_set_of_ns(const Index_Set& new_ns, NS_Rows& set_of_ns) {
  if (remove_redundants(new_ns, set_of_ns))
    set_of_ns.push_back(new_ns);
}
void
add_to_minimal_set_of_ns(Index_Set&& new_ns, NS_Rows& set_of_ns) {
  if (remove_redundants(new_ns, set_of_ns))
    set_of_ns.push_back(std::move(new_ns));
}

// For each ns in NS+-, recompute its support.
// If the new support is a singleton, move it into SK.
// Note: SK_Rows in std::vector<T>, where T in { Con, Gen }.
// Comments and parameter names are written for T = Gen.
// When T = Con, then the Boolean flag `is_strict_ineq'
// actually stores `is_point'.
template <typename SK_Rows>
void
move_ns(NS_Rows& ns_rows, SK_Rows& sk_rows, bool is_strict_ineq,
        Ranges& ns_ranges, NS_Rows& set_of_ns,
        const Index_Set& iset_neg, const Index_Set& iset_eq,
        const Sat& sat_src, const Sat& sat_dst) {
  Index_Partition<Range::POS_NEG> ns_pos_neg(ns_ranges);
  Dims to_be_del;
  for (auto i : ns_pos_neg) {
    assert(ns_ranges[i] == Range::POS_NEG);
    Index_Set& supp = ns_rows[i];
    support_closure(supp, sat_src, sat_dst);
    supp_projection(supp, is_strict_ineq, iset_neg, iset_eq);

    // Check if it can be moved into SK.
    if (!is_strict_ineq && supp.size() == 1) {
      // Change its only supporter's type into point.
      const auto supporter_index = *(supp.begin());
      promote_singleton(sk_rows[supporter_index]);
      to_be_del.push_back(i);
      // Consider next ns in pos_neg
      continue;
    }
    assert(supp.size() > 1);

    // Add `supp' to minimal `set_of_ns' to check its redundancy.
    add_to_minimal_set_of_ns(supp, set_of_ns);
  }
  // Remove ns that has been moved into SK.
  erase_using_sorted_indices(ns_rows, to_be_del);
  erase_using_sorted_indices(ns_ranges, to_be_del);
}

// Optimization function.
// Remove from `search_range' any x in `F_g' such that
// |sat_src(F_g)| = |sat_src(`ns_face' U {x})|.
void
remove_equivalent_generators(const Index_Set& F_g,
                             dim_type g,
                             const Bits& sat_F,
                             const Sat& sat_src,
                             Index_Set& search_range) {
  const dim_type F_g_ones = sat_F.count_ones_in_union(sat_src[g]);
  // Search only in search_range intersect F_g
  auto lookup_set = search_range & F_g;
  for (auto x : lookup_set) {
    const dim_type F_x_ones = sat_F.count_ones_in_union(sat_src[x]);
    if (F_g_ones == F_x_ones)
      search_range.reset(x);
  }
}

template <typename SK_Rows>
void
create_ns_with_enumerate_faces(const Index_Set& ns_face,
                               const Range sk_range,
                               const SK_Rows& sk_rows,
                               NS_Rows& set_of_ns,
                               const Index_Set& iset_pos,
                               const Index_Set& iset_eq,
                               const Index_Set& iset_neg,
                               const Sat& sat_src,
                               const Sat& sat_dst,
                               const bool is_strict_ineq) {
  // Search_range contains all the generators that can
  // introduce, with `ns_face', a face cut by the constraint.
  // Partition with Index_Set is used to make easier all
  // the supp_projections that will be needed in the end.
  // Besides, in case optimization is introduced,
  // using Index_Set for `search_range' should make easier removing
  // redundant generators from it, without editing `sk_ranges'

  assert(sk_range == Range::NEG || sk_range == Range::POS);
  Index_Set search_range
    = (sk_range == Range::NEG) ? iset_neg : iset_pos;
  const Bits sat_F = sat_all(ns_face, sat_src);

  for (auto g : search_range) {
    // Note: points can't be excluded for non-strict inequality.
    // Some faces shall be introduced by points that
    // will be lost during the projection onto NS=.
    if (is_strict_ineq && is_strict_ineq_or_point(sk_rows[g].type()))
      continue;

    // F_g = supp_g(`ns_face' U {g})
    Index_Set F_g(ns_face);
    F_g.set(g);
    support_closure(F_g, sat_src, sat_dst);

    // Although we are running on `search_range', removing elements
    // from it doesn't affect the end control.
    // FIXME: just for clarity, though, we could use a boolean mask
    // and keep `search_range' const (Index_Partition could be used too).
    remove_equivalent_generators(F_g, g, sat_F, sat_src, search_range);

    supp_projection(F_g, is_strict_ineq, iset_neg, iset_eq);

    // Redundancy check in set_of_ns
    add_to_minimal_set_of_ns(std::move(F_g), set_of_ns);
  } // end loop on g in search_range
}

template <typename SK_Rows>
void
create_ns_from_sk_with_adj_test(const SK_Rows& sk_rows,
                                const Ranges& sk_ranges,
                                NS_Rows& set_of_ns,
                                const Sat& sat_src,
                                const dim_type min_sat,
                                const dim_type max_sat) {
  const dim_type sk_rows_size = sk_rows.size();

  // Combine each point in SK= with *adjacent* elements in SK+.
  Index_Partition<Range::EQ> sk_eq(sk_ranges);
  Index_Partition<Range::POS> sk_pos(sk_ranges);
  for (auto eq : sk_eq) {
    if (!is_strict_ineq_or_point(sk_rows[eq].type()))
      continue;
    for (auto pos : sk_pos) {
      if (is_strict_ineq_or_point(sk_rows[pos].type()))
        continue;
      // `eq' is a point, `pos' is not a point.
      Bits new_satrow = Bits::get_union(sat_src[eq], sat_src[pos]);
      const dim_type new_satrow_ones = new_satrow.count_ones();
      if (!adj_test(0, sk_rows_size,
                    eq, pos, new_satrow, new_satrow_ones,
                    min_sat, max_sat, sat_src))
        // Adjacency test failed: consider next `pos'.
        continue;
      // Add {eq, pos} in new_ns, if non-redundant.
      Index_Set new_ns(std::make_pair(eq, pos));
      add_to_minimal_set_of_ns(std::move(new_ns), set_of_ns);
    }
  }
}

template<typename SK_Rows>
void
new_ns_from_set(const SK_Rows& sk_rows,
                NS_Rows& ns_dest,
                NS_Rows& set_of_ns,
                const Ranges& sk_ranges,
                const Ranges& ns_ranges,
                bool is_strict_ineq) {
  // For each face in `set_of_ns' add a new ns
  for (auto& new_ns : set_of_ns) {
    // Redundancy check with other points.

    // * Note 1: the "SP check" is much more effective if
    // done *after* the promote_singletons() call.
    // * Note 2: the "NS check" should be faster if done here
    // (instead of during the construction of `set_of_ns')
    // because only the projected support is checked.
    // Check if it might lose some redundant cases.

    bool to_be_inserted = true;

    // SP check
    // Check if in its updated support there is a point.
    for (auto new_ns_g : new_ns) {
      if (!is_strict_ineq_or_point(sk_rows[new_ns_g].type())) {
        continue;
      }
      // `new_ns_g' is a point:
      // If it is in Q_EQ with strict, it will become a closure point
      if (is_strict_ineq && sk_ranges[new_ns_g] == Range::EQ) {
        continue;
      }
      to_be_inserted = false;
      break;
    } // end loop on `new_ns_g'.
    if (!to_be_inserted) {
      // Found some SP in its support.
      // Consider next `new_ns'.
      continue;
    }

    // NS check
    // Check if in its updated support there is another NS.
    if (is_strict_ineq) {
      Index_Partition<Range::POS> ns_pos(ns_ranges);
      for (auto i : ns_pos) {
        if (subset_eq(ns_dest[i], new_ns)) {
          to_be_inserted = false;
          break;
        }
      }
    } else {
      Index_Partition<Range::EQ> ns_eq(ns_ranges);
      for (auto i : ns_eq) {
        if (subset_eq(ns_dest[i], new_ns)) {
          to_be_inserted = false;
          break;
        }
      }
    }

    if (to_be_inserted)
      ns_dest.push_back(new_ns);
  } // end loop on `new_ns' in `set_of_ns'.
}

// Performs the following combinations:
//  * points in Q- (SK- and NS-) are combined with SK+
//  if (strict_ineq)
//    * points in Q= (SK= and NS=) are combined with SK-
//      (maybe using adj test when processing SK=)
//  else
//    * points in Q+ (SK+ and NS+) are combined with SK-
//
// Combinations are added to `set_of_ns'; other parameters are all const.
template <typename SK_Rows>
void
create_ns(const SK_Rows& sk_rows, const NS_Rows& ns_rows,
          const Ranges& sk_ranges, const Ranges& ns_ranges,
          NS_Rows& set_of_ns,
          const Index_Set& iset_pos,
          const Index_Set& iset_eq,
          const Index_Set& iset_neg,
          const Sat& sat_src, const Sat& sat_dst,
          const bool is_strict_ineq,
          const dim_type min_sat, const dim_type max_sat) {
  // Combine points in Q- with elements in SK+.
  create_ns_from_Q_with_face_enum<Range::NEG>(Range::POS,
                                              sk_rows, ns_rows,
                                              sk_ranges, ns_ranges,
                                              set_of_ns,
                                              iset_pos, iset_eq, iset_neg,
                                              sat_src, sat_dst,
                                              is_strict_ineq);
  if (is_strict_ineq) {
    // Combine points in Q= with elements in SK+.
    create_ns_from_eq_points(sk_rows, ns_rows, sk_ranges, ns_ranges,
                             set_of_ns, iset_pos, iset_eq, iset_neg,
                             sat_src, sat_dst, min_sat, max_sat);
  } else {
    // Combine points in Q+ with elements in SK-.
    create_ns_from_Q_with_face_enum<Range::POS>(Range::NEG,
                                                sk_rows, ns_rows,
                                                sk_ranges, ns_ranges,
                                                set_of_ns,
                                                iset_pos, iset_eq, iset_neg,
                                                sat_src, sat_dst, false);
  }
}

// Process points saturating strict ineq
// (resp., strict inequalities saturated by point).
template <typename SK_Rows>
void
process_points_sat_strict(SK_Rows& sk_rows, NS_Rows& ns_rows,
                          const Ranges& sk_ranges, const Ranges& ns_ranges,
                          const Sat& sat_src, Sat& sat_dst,
                          const dim_type min_sat, const dim_type max_sat) {
  assert(sk_ranges.size() == sk_rows.size());
  assert(ns_ranges.size() == ns_rows.size());
  // Update sat_dst.
  sat_dst = sat_src.transpose();

  Index_Set iset_neg, iset_eq, iset_pos;
  sk_ranges_to_index_sets(sk_ranges, iset_neg, iset_eq, iset_pos);

  NS_Rows set_of_ns;
  // Combine points in Q= (SK= and NS=) with NS+.
  // All new points are in NS+, so that NS= is not affected
  // (this is important, as it will be used just before returning).
  // Since `ns_ranges' will not be used after this and until
  // the end of the main iteration, no need to update it.
  create_ns_from_eq_points(sk_rows, ns_rows, sk_ranges, ns_ranges,
                           set_of_ns, iset_pos, iset_eq, iset_neg,
                           sat_src, sat_dst, min_sat, max_sat);
  // Actually add new non-redundant NP.
  new_ns_from_set(sk_rows, ns_rows, set_of_ns, sk_ranges, ns_ranges, true);

  // Each point in SK= becomes a closure point.
  Index_Partition<Range::EQ> sk_eq(sk_ranges);
  points_become_closure_points(sk_eq, sk_rows);
  // Each point in NS= is removed.
  Index_Partition<Range::EQ> ns_eq(ns_ranges);
  erase_using_sorted_indices(ns_rows, ns_eq);
}

template <typename SK_Rows>
void
create_ns_with_ex_sing(SK_Rows& sk_rows,
                       NS_Rows& ns_rows,
                       const Ranges& sk_ranges,
                       const Ranges& ns_ranges,
                       dim_type ex) {
  assert(sk_ranges.size() == sk_rows.size());
  assert(ns_ranges.size() == ns_rows.size());

  // For each point in SK=, create a new point in NS+ with ex.
  Index_Partition<Range::EQ> sk_eq(sk_ranges);
  for (auto eq : sk_eq) {
    if (!is_strict_ineq_or_point(sk_rows[eq].type()))
      continue;
    Index_Set new_ns(std::make_pair(eq, ex));
    ns_rows.push_back(new_ns);
  }
  // For each NP in NS=, add ex to its supporters
  Index_Partition<Range::EQ> ns_eq(ns_ranges);
  for (auto eq : ns_eq)
    ns_rows[eq].set(ex);

  // Each point in SK= becomes a closure point.
  points_become_closure_points(sk_eq, sk_rows);
  // Each element in NS= has been moved in NS+ adding `ex',
  // hence there are no more elements in NS= to be removed.
}

/* Note: used for debugging only; no need to optimize. */
template <typename SK_Rows>
bool
check_redundant_ns_by_sk(const NS_Rows& ns_rows,
                         const SK_Rows& sk_rows) {
  Index_Set skel_ns;
  for (auto i : bwd_index_range(sk_rows)) {
    if (is_strict_ineq_or_point(sk_rows[i].type()))
      skel_ns.set(i);
  }
  for (const auto& ns : ns_rows) {
    if (ns.intersects(skel_ns))
      return false;
  }
  return true;
}

/* Note: used for debugging only; no need to optimize. */
bool
check_redundant_ns_by_ns(const NS_Rows& ns_rows) {
  for (auto i : bwd_index_range(ns_rows)) {
    for (auto j : bwd_index_range(ns_rows)) {
      if (i == j) continue;
      if (subset_eq(ns_rows[i], ns_rows[j]))
        return false;
    }
  }
  return true;
}

/* Note: used for debugging only; no need to optimize. */
bool
check_ns_rows(const NS_Rows& ns_rows) {
  return all_of(ns_rows, [](const Index_Set& ns) { return ns.size() > 1; });
}

template <typename Part>
inline bool
has_strict_in_range(const Part& part, const Cons& rows) {
  return any_of(part, [&rows](dim_type i) {
                        return rows[i].is_strict_inequality();
                      });
}

template <typename Src_Row, typename Dst_Sys>
void
process_violating_singular(const dim_type ex,
                           const Src_Row& src,
                           Dst_Sys& dst,
                           Sat& sat_src,
                           Ranges& sk_ranges,
                           Ranges& ns_ranges) {
  auto& sg_dst = dst.sing_rows;
  auto& sk_dst = dst.sk_rows;
  auto& ns_dst = dst.ns_rows;
  // Map the violating line into a ray.
  auto& sg_ex = sg_dst[ex];
  Integer sp_ex;
  sp::assign(sp_ex, src, sg_ex);
  make_non_singular(sg_ex, sp_ex);
  // Map all other geometric rows into SK=.
  // Note: we do not create/move rows, hence ns_rows is not affected.
  combine_with_ex_sing(ex, sp_ex, src, sg_dst, sk_dst);

  if (is_singular(src)) {
    // Just remove it (do not update sat_src).
    sg_dst.erase(sg_dst.begin() + ex);
    return;
  }

  // Since sg_ex is no longer singular, move it to sk_dst.
  const auto new_ex = num_rows(sk_dst);
  sk_dst.push_back(std::move(sg_ex));
  sg_dst.erase(sg_dst.begin() + ex);

  // Update `sat': first add a new column ...
  assert(!is_singular(src));
  sat_src.add_cols(1);
  // ... then add a new row for the ex singular moved into sk_dst;
  // the ex singular saturates all but the newly added column.
  Bits sat_ex;
  sat_ex.set(sat_src.num_cols() - 1);
  sat_src.add_row(std::move(sat_ex));

  if (is_strict_ineq_or_point(src.type())) {
    // Create partitions for sk_dst and ns_dst.
    sk_ranges.assign(sk_dst.size(), Range::EQ);
    sk_ranges[new_ex] = Range::POS;
    ns_ranges.assign(ns_dst.size(), Range::EQ);
    create_ns_with_ex_sing(sk_dst, ns_dst, sk_ranges, ns_ranges, new_ex);
  }
}

inline bool
maybe_has_positivity(const Cons& cons) {
  return any_of(cons, std::mem_fn(&Con::is_tautological));
}
inline bool
maybe_has_positivity(const Gens&) {
  return false;
}
inline bool
maybe_is_ray(const Con&) {
  return false;
}
inline bool
maybe_is_ray(const Gen& g) {
  return g.is_ray();
}
inline bool
maybe_has_ray(const Cons&) {
  return false;
}
inline bool
maybe_has_ray(const Gens& gens) {
  return any_of(gens, std::mem_fn(&Gen::is_ray));
}

inline bool
can_use_gen_to_con_fast_path(const Con&, const Gens&, const NS_Rows&,
                             bool, bool&, bool&) {
  /* The fast path is only used for gen-to-con conversions */
  return false;
}
inline bool
can_use_gen_to_con_fast_path(const Gen& src,
                             const Cons& sk_dst,
                             const NS_Rows& ns_dst,
                             bool /*has_rays*/,
                             bool& has_positivity,
                             bool& has_maximal_efc) {
  FAST_PATH_STATS_INCR_CALLS;
  // Fast path disabled when processing a closure point.
  if (src.is_closure_point())
    return false;
  // Fast path disabled when we have a (non-trivial) strict inequality.
  const auto num_ns = ns_dst.size();
  if (num_ns > 1)
    return false;
  if (any_of(sk_dst,
             [](const Con& c) {
               return c.is_strict_inequality() && not c.is_tautological();
             }))
    return false;

  // Check for (skel) positivity constraint.
  has_positivity = any_of(sk_dst, std::mem_fn(&Con::is_tautological));
  if (has_positivity) {
    has_maximal_efc = false;
    if (num_ns == 0) {
      FAST_PATH_STATS_INCR_TRUE;
      return true;
    } else
      return false;
  }
  // Check for (non-skel) maximal efc.
  has_maximal_efc = (num_ns == 1) && (ns_dst[0].size() == num_rows(sk_dst));
  if (has_maximal_efc || (num_ns == 0)) {
    FAST_PATH_STATS_INCR_TRUE;
    return true;
  }

  assert(not has_positivity && not has_maximal_efc && (num_ns == 1));
  // The (one and only) support in ns_dst could be the rfc,
  // however, we miss a way to cheaply check if this is the case;
  // hence, we conservatively avoid the fast path.
  FAST_PATH_STATS_INCR_DOUBT;
  return false;
  /*
    Note: if we could be sure that ns_dst[0] is the rfc,
    then we could do the following:
      return not src.is_ray();
    because we need to update the rfc only when adding a ray.
  */
}

inline bool
should_call_create_ns(const Con& src, const Gens& sk_dst,
                      const Ranges& sk_ranges, const Ranges&) {
  Index_Partition<Range::EQ> sk_eq(sk_ranges);
  return src.is_strict_inequality() || can_have_nonskel(sk_eq, sk_dst);
}

inline bool
should_call_create_ns(const Gen& src, const Cons& sk_dst,
                      const Ranges& sk_ranges, const Ranges& ns_ranges) {
  const bool is_point = src.is_point();
  Index_Partition<Range::EQ> sk_eq(sk_ranges);
  if (!is_point && !can_have_nonskel(sk_eq, sk_dst))
    return false;
  Index_Partition<Range::NEG> sk_neg(sk_ranges);
  Index_Partition<Range::NEG> ns_neg(ns_ranges);
  if (!ns_neg.empty() || has_strict_in_range(sk_neg, sk_dst))
    return true;
  Index_Partition<Range::EQ> ns_eq(ns_ranges);
  if (is_point && (!ns_eq.empty() || has_strict_in_range(sk_eq, sk_dst)))
    return true;
  Index_Partition<Range::POS> sk_pos(sk_ranges);
  Index_Partition<Range::POS> ns_pos(ns_ranges);
  return !is_point
    && (!ns_pos.empty() || has_strict_in_range(sk_pos, sk_dst));
}

/* Only for gen-to-con conversion */
void
ensure_strict_pos(Cons& sk_dst, NS_Rows& ns_dst) {
  // At most one tautological constraint in sk_dst.
  assert(std::count_if(sk_dst.begin(), sk_dst.end(),
                       std::mem_fn(&Con::is_tautological)) <= 1);
  for (auto pos_i : bwd_index_range(sk_dst)) {
    auto& pos = sk_dst[pos_i];
    if (!pos.is_tautological())
      continue;
    // Found a tautological constraint (it has to be the only one).
    if (pos.is_strict_inequality())
      return;
    // Make it strict.
    pos.set_type(Con::STRICT_INEQUALITY);
    // At most one support in ns_dst contains index `pos_i'.
    assert(std::count_if(ns_dst.begin(), ns_dst.end(),
                         [pos_i](const Index_Set& ns) {
                           return ns.test(pos_i);
                         }) <= 1);
    // Find if there is an efc holding the positivity in ns_dst.
    auto efc_i = std::find_if(ns_dst.begin(), ns_dst.end(),
                              [pos_i](const Index_Set& ns) {
                                return ns.test(pos_i);
                              });
    // If found, remove it.
    if (efc_i != ns_dst.end())
      ns_dst.erase(efc_i);
    return;
  }
}

/* Note: used for debugging only. */
/* Returns false if the first positivity constraint found is non-strict.
   Note: we only check the first, as sk_dst is assumed non-redundant. */
bool
check_strict_pos(const Cons& sk_dst) {
  auto pos_i = std::find_if(sk_dst.begin(), sk_dst.end(),
                            std::mem_fn(&Con::is_tautological));
  return pos_i == sk_dst.end() || pos_i->is_strict_inequality();
}

/**
 * Make sure that there is a cutter and it is not redundant.
 * Not checking if there is the correct rfc,
 * as this control would need sat_c too.
 */
bool
check_cutters(const NS_Rows& ns_dst,
              const Cons& sk_dst,
              const Gens& sk_src,
              const Sat& sat_g) {
  assert(sat_g.num_cols() == num_rows(sk_src));

  // Check for at least one face cutter (skel or non-skel).
  bool has_sk_cutter = has_strict_ineq(sk_dst);
  bool has_ns_cutter = !ns_dst.empty();

  if (!has_sk_cutter && !has_ns_cutter)
    return false;

  // The positivity constraint, if present, must be strict.
  auto sk_pos = std::find_if(sk_dst.cbegin(), sk_dst.cend(),
                             std::mem_fn(&Con::is_tautological));
  bool has_sk_pos = (sk_pos != sk_dst.cend());
  if (has_sk_pos && !sk_pos->is_strict_inequality())
    return false;

  // Positivity and efc are mutually exclusive.
  auto efc = std::find_if(ns_dst.cbegin(), ns_dst.cend(),
                          [&sat_g, &sk_src](const Index_Set& ns) {
                            return is_empty_face_cutter(ns, sat_g, sk_src);
                          });
  bool has_efc = (efc != ns_dst.cend());
  if (has_efc && has_sk_pos)
    return false;

  // In principle, here we should only check for efc non-redundancy.
  // We perform a complete non-skel non-redundancy test
  // (this code is only executed in debugging mode).
  return check_redundant_ns_by_sk(ns_dst, sk_dst)
    && check_redundant_ns_by_ns(ns_dst);
}

template <bool con_to_gen, typename Src, typename Dst_Sys>
Conv_Iter_Result
conv_iter(dim_type space_dim, dim_type src_num_sing,
          bool has_rays,
          const Src& src, Dst_Sys& dst,
          Sat& sat_src, Sat& sat_dst,
          Integers& scal_prods, NS_Rows& set_of_ns,
          Ranges& sk_ranges, Ranges& ns_ranges,
          bool is_support) {
  auto& sg_dst = dst.sing_rows;
  auto& sk_dst = dst.sk_rows;
  auto& ns_dst = dst.ns_rows;

  // Check redundant non-skel generators.
  assert(check_redundant_ns_by_sk(ns_dst, sk_dst));
  assert(check_redundant_ns_by_ns(ns_dst));

  // Check if all singular rows saturate src.
  const dim_type ex = check_sing_rows(src, sg_dst);
  if (ex < num_rows(sg_dst)) {
    process_violating_singular(ex, src, dst, sat_src,
                               sk_ranges, ns_ranges);
    return Conv_Iter_Result::OK;
  }

  // Here all lines saturate the constraint.
  sp::compute_prods(src, sk_dst, scal_prods);

  // Create partition for skeleton rows.
  sk_partition(scal_prods, sk_ranges);
  Index_Partition<Range::NEG> sk_neg(sk_ranges);
  Index_Partition<Range::EQ> sk_eq(sk_ranges);
  Index_Partition<Range::POS> sk_pos(sk_ranges);

  const bool singular = is_singular(src);
  const bool strict_or_point = is_strict_ineq_or_point(src.type());

  // Check if SK- is empty.
  if (sk_neg.empty()) {

    if (singular) {
      if (sk_pos.empty())
        return Conv_Iter_Result::REDUNDANT;
      if (con_to_gen && sk_eq.empty())
        return Conv_Iter_Result::EMPTY;
      // Singular with SK- empty and SK= non-empty: remove SK+ and NS+.
      // Delete rows in NS+.
      sk_to_ns_ranges(sk_ranges, ns_dst, ns_ranges);
      Index_Partition<Range::POS> ns_pos(ns_ranges);
      erase_using_sorted_indices(ns_dst, ns_pos);
      // Delete rows in SK+.
      remove_rows(sk_pos, sk_dst, ns_dst, sat_src);
      return Conv_Iter_Result::OK;
    }

    if (!strict_or_point) {
      if (sk_eq.empty() || !is_support)
        // src is redundant and makes redundant every ns including it.
        return Conv_Iter_Result::REDUNDANT;
      else {
        // Update sat matrix.
        assert(!is_singular(src));
        sat_src.add_cols(1);
        update_last_sat_column(sat_src, sk_pos);
        // ENEA: FIXME: consider if another return value may be useful.
        // This skel element is redundant, but it can not be immediately
        // removed because it may appear in pending.ns_rows.
        return Conv_Iter_Result::OK;
      }
    }

    assert(strict_or_point);
    if (con_to_gen && sk_pos.empty())
      return Conv_Iter_Result::EMPTY;
    if (sk_eq.empty())
      return Conv_Iter_Result::REDUNDANT;
    sk_to_ns_ranges(sk_ranges, ns_dst, ns_ranges);
    // Update sat.
    assert(!is_singular(src));
    sat_src.add_cols(1);
    update_last_sat_column(sat_src, sk_pos);
    const dim_type min_sat = space_dim - sg_dst.size() - 1;
    const dim_type max_sat = sat_src.num_cols() + src_num_sing;
    process_points_sat_strict(sk_dst, ns_dst, sk_ranges, ns_ranges,
                              sat_src, sat_dst, min_sat, max_sat);

    return Conv_Iter_Result::OK;
  } // end of case analysis for SK- empty

  // Here SK- is not empty: check if SK+ is empty.
  if (sk_pos.empty()) {
    if (con_to_gen && strict_or_point)
      return Conv_Iter_Result::EMPTY;
    if (con_to_gen && sk_eq.empty())
      return Conv_Iter_Result::EMPTY;

    // SK+ empty, SK= non-empty: remove rows in SK- and NS-.
    // Delete rows in NS-.
    sk_to_ns_ranges(sk_ranges, ns_dst, ns_ranges);
    Index_Partition<Range::NEG> ns_neg(ns_ranges);
    erase_using_sorted_indices(ns_dst, ns_neg);
    // Delete rows in SK-.
    remove_rows(sk_neg, sk_dst, ns_dst, sat_src);
    if (!singular) {
      // Update `sat_src': new column all made of zeroes.
      // FIXME: we have identified an implicit singular src row.
      // Consider returning a different result value to signal
      // the caller that it can be moved into src.sing_rows
      assert(!is_singular(src));
      sat_src.add_cols(1);
    }
    return Conv_Iter_Result::OK;
  } // end of case analysis for SK- non-empty and SK+ empty.

  // Here SK+ and SK- are both non-empty.
  assert(!sk_pos.empty() && !sk_neg.empty());

  bool has_positivity = false;
  bool has_maximal_efc = false;
  bool use_fast_path
    = can_use_gen_to_con_fast_path(src, sk_dst, ns_dst,
                                   has_rays,
                                   has_positivity,
                                   has_maximal_efc);

  // Combine SK+ and SK- *adjacent* generators.
  { // NOTE: max_sat has to be computed *before* updating sat,
    // since otherwise the quick_non_adj test is not correct.
    const dim_type min_sat = space_dim - sg_dst.size() - 1;
    const dim_type max_sat = sat_src.num_cols() + src_num_sing;
    combine_sk_rows(sk_ranges, sk_dst, scal_prods, sat_src,
                    strict_or_point, min_sat, max_sat);
  }

  if (!singular) {
    // Update sat.
    assert(!is_singular(src));
    sat_src.add_cols(1);
    update_last_sat_column(sat_src, sk_pos);
  }

  // Create partition for non-skeleton rows.
  // Note: computing `iset_eq' *after* combine_sk_rows()
  // means it is up-to-date for later calls to move_ns().
  Index_Set iset_neg, iset_eq, iset_pos;
  sk_ranges_to_index_sets(sk_ranges, iset_neg, iset_eq, iset_pos);
  ns_partition(ns_dst, iset_neg, iset_pos, ns_ranges);

  // Refine the test for the gen-to-con fast path.
  if (use_fast_path && not ns_dst.empty() && not has_maximal_efc) {
    assert(num_rows(ns_dst) == 1);
    // ns_dst[0] is a (non maximal efc) strict inequality;
    // if it is in POS_NEG, we cannot take the fast path;
    if (ns_ranges[0] == Range::POS_NEG)
      use_fast_path = false;
  }

  if (use_fast_path) {
    FAST_PATH_STATS_INCR_EXEC_FAST;
    // gen-to-con fast path for the cases when
    // there is no need to transpose the saturation matrix.
    assert(ns_dst.size() <= 1);
    // If there is the positivity then there is nothing to do.
    if (!has_positivity && !maybe_has_positivity(sk_dst)) {
      ns_ranges[0] = singular ? Range::EQ : Range::POS;
      if (is_strict_ineq_or_point(src.type())) {
        // src is point.
        if (has_maximal_efc) {
          // Just keep the efc maximal.
          auto& ns_efc = ns_dst[0];
          ns_efc.set_until(sk_dst.size());
        }
        // If it was the rfc, nothing to do.
      } else {
        // src is ray or line.
        if (has_maximal_efc) {
          // It was closed, and it is the first ray / line.
          // Transform it in rfc:
          ns_dst[0] = iset_eq;
        } else {
          assert(is_singular(src));
          assert(ns_dst.size() == 1);
          // Nothing to do.
        }
      }
    }
  } else {
    FAST_PATH_STATS_INCR_EXEC_SLOW;
    // Dealing with non-skel elements (not the fast path).
    // New non-skel elements will be placed here.
    set_of_ns.clear();
    bool sat_dst_is_up_to_date = false;

    // Check if we need to call `move_ns'.
    Index_Partition<Range::POS_NEG> ns_pos_neg(ns_ranges);
    if (!ns_pos_neg.empty()) {
      assert(!sat_dst_is_up_to_date);
      sat_dst = sat_src.transpose();
      sat_dst_is_up_to_date = true;
      move_ns(ns_dst, sk_dst, strict_or_point,
              ns_ranges, set_of_ns, iset_neg, iset_eq,
              sat_src, sat_dst);
    }

    // Check if we need to call `create_ns'.
    const bool calling_create_ns
      = should_call_create_ns(src, sk_dst, sk_ranges, ns_ranges);
    if (calling_create_ns) {
      if (!sat_dst_is_up_to_date) {
        sat_dst = sat_src.transpose();
        sat_dst_is_up_to_date = true;
      }
      // Combine ns with *non-adjacent* rows.
      const dim_type min_sat = space_dim - sg_dst.size() - 1;
      const dim_type max_sat = sat_src.num_cols() + src_num_sing;
      create_ns(sk_dst, ns_dst, sk_ranges, ns_ranges,
                set_of_ns, iset_pos, iset_eq, iset_neg,
                sat_src, sat_dst, strict_or_point,
                min_sat, max_sat);
    }

    // If we updated `sat_dst', then we called `move_ns' or `create_ns',
    // hence we need to add `set_of_ns' to `ns_dest'.
    if (sat_dst_is_up_to_date) {
      // Actually add new non-redundant NP.
      new_ns_from_set(sk_dst, ns_dst, set_of_ns, sk_ranges, ns_ranges,
                      strict_or_point);

      if (strict_or_point)
        points_become_closure_points(sk_eq, sk_dst);
      else {
        // Note: if the assertion below ever turns out to be false
        // due to the presence of a singleton ns,
        // then we will have to call helper function
        //   promote_singletons(ns_dst, sk_dst);
        // and *also* erase the promoted ns from ns_ranges.
        assert(check_ns_rows(ns_dst));
      }
    }
  }

  if (is_singular(src)) {
    // Modify sk_ranges and ns_ranges so that
    // rows in SK+ and NS+ will get removed.
    std::replace(sk_ranges.begin(), sk_ranges.end(),
                 Range::POS, Range::NEG);
    std::replace(ns_ranges.begin(), ns_ranges.end(),
                 Range::POS, Range::NEG);
  } else if (strict_or_point) {
    // Modify ns_ranges so that rows in NS= will get removed.
    std::replace(ns_ranges.begin(), ns_ranges.end(),
                 Range::EQ, Range::NEG);
  }
  // Modify ns_ranges so that rows in NS+- will get removed.
  std::replace(ns_ranges.begin(), ns_ranges.end(),
               Range::POS_NEG, Range::NEG);

  // Delete rows in NS-.
  Index_Partition<Range::NEG> ns_neg(ns_ranges);
  erase_using_sorted_indices(ns_dst, ns_neg);
  // Delete rows in SK-.
  remove_rows(sk_neg, sk_dst, ns_dst, sat_src);

  return Conv_Iter_Result::OK;
}

template <typename Row>
void
move_to_end(std::vector<Row>& src, std::vector<Row>& dst) {
  dst.reserve(src.size() + dst.size());
  for (auto& row : src)
    dst.push_back(std::move(row));
  src.clear();
}

/* Note: these checks are just an optimization */
// If the efc is not skipped here, simplify will remove it.
// Consider if it is worth it.
bool
is_efc_from_sat(const Index_Set&, const Cons&) {
  /* Only for con2gen conversion. */
  PPLITE_UNREACH;
  return false;
}
inline bool
is_efc_from_sat(const Index_Set& eq,
                const Gens& sk) {
  return all_of(eq, [&sk](const dim_type i){return sk[i].is_ray();});
}
template <typename Dst>
bool
is_redundant(bool con_to_gen,
             const Index_Set& eq,
             const Dst& dst) {
  if (eq.empty())
    return true;
  if (!con_to_gen)
    return false;
  // In con_to_gen conversion exclude the efc too.
  // Since it comes from another polyhedron, it
  // may not be in the maximal representation.
  return is_efc_from_sat(eq, dst);
}

template <bool con_to_gen, typename Src_Sys, typename Dst_Sys>
Conv_Iter_Result
conv_sing(dim_type space_dim,
          Src_Sys& src, Src_Sys& pending, Dst_Sys& dst,
          bool& has_rays,
          Sat& sat_src, Sat& sat_dst,
          Integers& scal_prods, NS_Rows& set_of_ns,
          Ranges& sk_ranges, Ranges& ns_ranges) {
  for (auto k : index_range(pending.sing_rows)) {
    auto& src_k = pending.sing_rows[k];
    OUTPUT_COUNTERS_IF_NOISY("sing", k, dst);
    assert(maybe_check_cutters(dst.ns_rows, dst.sk_rows,
                               src.sk_rows, sat_src));
    auto res = conv_iter<con_to_gen>(space_dim, num_rows(src.sing_rows),
                                     has_rays,
                                     src_k, dst, sat_src, sat_dst,
                                     scal_prods, set_of_ns,
                                     sk_ranges, ns_ranges,
                                     false);
    switch (res) {
    case Conv_Iter_Result::OK:
      src.sing_rows.push_back(std::move(src_k));
      // Normalize positivity (if present) to be a strict constraint.
      maybe_ensure_strict_pos(dst.sk_rows, dst.ns_rows);
      break;
    case Conv_Iter_Result::EMPTY:
      return res;
    case Conv_Iter_Result::REDUNDANT:
      // Nothing to do.
      break;
    }
  }
  // pending.sing_rows now contains moved/redundant rows only.
  pending.sing_rows.clear();
  return Conv_Iter_Result::OK;
}

template <bool con_to_gen, typename Src_Sys, typename Dst_Sys>
Conv_Iter_Result
conv_skel(dim_type space_dim,
          Src_Sys& src, Src_Sys& pending, Dst_Sys& dst,
          bool& has_rays,
          Sat& sat_src, Sat& sat_dst,
          Integers& scal_prods, NS_Rows& set_of_ns,
          Ranges& sk_ranges, Ranges& ns_ranges) {
#if PPLITE_SORT_PENDING
  Dims sorted;
  { // Populate sorted with indices on pending.sk_rows
    // and then (indirectly) sort them according to lexicographic order.
    // Clock clock;
    const auto& sk_rows = pending.sk_rows;
    auto cmp = [&sk_rows](dim_type i, dim_type j) {
      return compare(sk_rows[i], sk_rows[j]) < 0;
    };
    sorted.resize(sk_rows.size());
    std::iota(sorted.begin(), sorted.end(), 0);
    std::sort(sorted.begin(), sorted.end(), cmp);
    // std::cerr << "\nsorting time: ";
    // clock.print_elapsed(std::cerr);
  }
#endif // PPLITE_SORT_PENDING
  Index_Set redundant;
  assert(pending.sing_rows.empty());
  // Store old src skel size, to be used after the loop.
  const dim_type src_num_sk = src.sk_rows.size();
  const dim_type src_num_sing = src.sing_rows.size();
  for (auto k : index_range(pending.sk_rows)) {
    const bool is_supp = any_of(pending.ns_rows,
				[k](const Index_Set& ns) {
                                  return ns.test(k);
				});
    // Note: using a non-const reference, as we will std::move() it.
#if PPLITE_SORT_PENDING
    auto& src_k = pending.sk_rows[sorted[k]];
#else
    auto& src_k = pending.sk_rows[k];
#endif
    OUTPUT_COUNTERS_IF_NOISY("skel", k, dst);
    assert(maybe_check_cutters(dst.ns_rows, dst.sk_rows,
                               src.sk_rows, sat_src));
    auto res = conv_iter<con_to_gen>(space_dim, src_num_sing,
                                     has_rays,
                                     src_k, dst, sat_src, sat_dst,
                                     scal_prods, set_of_ns,
                                     sk_ranges, ns_ranges,
                                     is_supp);
    switch (res) {
    case Conv_Iter_Result::OK:
      if (!has_rays && maybe_is_ray(src_k))
        has_rays = true;
      src.sk_rows.push_back(std::move(src_k));
      // Normalize positivity (if present) to be a strict constraint.
      maybe_ensure_strict_pos(dst.sk_rows, dst.ns_rows);
      break;
    case Conv_Iter_Result::EMPTY:
      return res;
    case Conv_Iter_Result::REDUNDANT:
      redundant.set(k);
      break;
    }
  }
  // pending.sk_rows now contains moved/redundant rows only.
  pending.sk_rows.clear();

#if PPLITE_SORT_PENDING
  // Recompute redundant according to `sorted'.
  // Index_Set inv_red;
  // for (auto i : redundant)
  //   inv_red.set(sorted[i]);
  // redundant = std::move(inv_red);
  // Set up a lambda for remapping supports according to sorted.
  auto sz = num_rows(sorted);
  Dims inv(sz, not_a_dim());
  for (auto i : range(sz)) {
    assert(inv[sorted[i]] == not_a_dim());
    inv[sorted[i]] = i;
  }
  auto remap = [&inv](const Index_Set& in) {
    Index_Set out;
    for (auto i : in) out.set(inv[i]);
    return out;
  };
#endif // PPLITE_SORT_PENDING

  // Adjust indices in pending.ns_rows.
  Index_Set ns_redundant;
  for (auto i : bwd_index_range(pending.ns_rows)) {
    auto& ns = pending.ns_rows[i];
#if PPLITE_SORT_PENDING
    ns = remap(ns);
#endif
    if (ns.intersects(redundant))
      ns_redundant.set(i);
    else {
      ns.remove_all(redundant);
      ns >>= src_num_sk;
    }
  }
  // Remove early detected redundant non-skel rows.
  if (!ns_redundant.empty())
    erase_using_sorted_indices(pending.ns_rows, ns_redundant);

  return Conv_Iter_Result::OK;
}

template <bool con_to_gen, typename Src_Sys, typename Dst_Sys>
void
conv_nonskel(dim_type space_dim,
             Src_Sys& src, Src_Sys& pending, Dst_Sys& dst,
             Sat& sat_src, Sat& sat_dst,
             Ranges& sk_ranges, Ranges& ns_ranges) {
  auto& sk_rows = dst.sk_rows;
  auto& ns_rows = dst.ns_rows;
  // FIXME: may be already up-to-date.
  sat_dst = sat_src.transpose();
  for (auto k : index_range(pending.ns_rows)) {
    auto& src_k = pending.ns_rows[k];
    Index_Set iset_pos { sat_all(src_k, sat_dst) };
    Index_Set iset_eq = iset_pos;
    iset_eq.complement_until(sat_dst.num_cols());
    if (is_redundant(con_to_gen, iset_eq, sk_rows))
      continue;
    Index_Set iset_neg;
    // Compute sk_ranges.
    sk_ranges.resize(sk_rows.size());
    for (auto i : bwd_index_range(sk_rows))
      sk_ranges[i] = iset_pos.test(i) ? Range::POS : Range::EQ;
    // Compute ns_ranges.
    ns_ranges.resize(ns_rows.size());
    for (auto i : bwd_index_range(ns_rows)) {
      bool has_pos = ns_rows[i].intersects(iset_pos);
      ns_ranges[i] = has_pos ? Range::POS : Range::EQ;
    }
    // This is the ending part of `process_point_sat_strict'.
    NS_Rows set_of_ns;
    assert(pending.sing_rows.empty());
    const dim_type min_sat = space_dim - dst.sing_rows.size() - 1;
    const dim_type max_sat = sat_src.num_cols() + num_rows(src.sing_rows);
    create_ns_from_eq_points(sk_rows, ns_rows, sk_ranges, ns_ranges,
                             set_of_ns, iset_pos, iset_eq, iset_neg,
                             sat_src, sat_dst, min_sat, max_sat);
    new_ns_from_set(sk_rows, ns_rows,
                    set_of_ns, sk_ranges, ns_ranges, true);
    Index_Partition<Range::EQ> sk_eq(sk_ranges);
    points_become_closure_points(sk_eq, sk_rows);
    Index_Partition<Range::EQ> ns_eq(ns_ranges);
    erase_using_sorted_indices(ns_rows, ns_eq);
    // Finished processing src_k: move it into src.
    src.ns_rows.push_back(std::move(src_k));
  }
  // pending.ns_rows now contains moved/redundant rows only.
  pending.ns_rows.clear();
  // CHECKME: No need to call maybe_ensure_strict_pos().
}

///////////////////////////////////////////////////////////////////////

// Function template: conversion.

// Template definition.
template <bool con_to_gen, typename Src_Sys, typename Dst_Sys>
bool
conversion(dim_type space_dim,
           Src_Sys& src, Src_Sys& pending,
           Dst_Sys& dst, Sat& sat_src) {
  // Allocate once for all.
  Integers scal_prods;
  Ranges sk_ranges, ns_ranges;
  Sat sat_dst;
  NS_Rows set_of_ns;
  bool has_rays = maybe_has_ray(src.sk_rows);

  FAST_PATH_STATS_RESET;

  // Process pending sing rows.
  if (!pending.sing_rows.empty()) {
    auto res = conv_sing<con_to_gen>(space_dim, src, pending, dst,
                                     has_rays,
                                     sat_src, sat_dst,
                                     scal_prods, set_of_ns,
                                     sk_ranges, ns_ranges);
    if (res == Conv_Iter_Result::EMPTY)
      return true;
  }

  // Process pending skel rows.
  if (!pending.sk_rows.empty()) {
    auto res = conv_skel<con_to_gen>(space_dim, src, pending, dst,
                                     has_rays,
                                     sat_src, sat_dst,
                                     scal_prods, set_of_ns,
                                     sk_ranges, ns_ranges);
    if (res == Conv_Iter_Result::EMPTY)
      return true;
  }

  // Process pending non-skel rows.
  if (!pending.ns_rows.empty()) {
    conv_nonskel<con_to_gen>(space_dim, src, pending, dst,
                             sat_src, sat_dst, sk_ranges, ns_ranges);
  }

  // Check for emptiness.
  if (con_to_gen && dst.ns_rows.empty() && !contains_a_point(dst.sk_rows))
    return true;

  FAST_PATH_STATS_PRINT;

  assert(pending.empty());
  assert(num_rows(src.sk_rows) == sat_src.num_cols());
  assert(num_rows(dst.sk_rows) == sat_src.num_rows());
  assert(maybe_check_cutters(dst.ns_rows, dst.sk_rows,
                             src.sk_rows, sat_src));
  return false;
}

// Explicit instantiations.
template
bool conversion<true, Con_Sys, Gen_Sys>(dim_type, Con_Sys&, Con_Sys&,
                                        Gen_Sys&, Sat&);
template
bool conversion<false, Gen_Sys, Con_Sys>(dim_type, Gen_Sys&, Gen_Sys&,
                                         Con_Sys&, Sat&);

} // namespace detail
} // namespace pplite
