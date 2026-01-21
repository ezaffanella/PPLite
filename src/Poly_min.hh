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

#ifndef pplite_Poly_min_hh
#define pplite_Poly_min_hh 1

#include "globals.hh"
#include "support_utils.hh"
#include "utils.hh"
#include "Bits.hh"
#include "Con.hh"
#include "Gen.hh"
#include "Index_Partition.hh"
#include "Integer.hh"
#include "Linear_Expr.hh"
#include "Poly.hh"
#include "Sat.hh"
#include "Scalar_Prod.hh"

#include <algorithm>
#include <cassert>
#include <iterator>

namespace pplite {

/* Extend namespace sp by adding methods for Poly_Impl::Sys<T> */
namespace sp {
bool
satisfied_by_all(const Topol topol, const Con& c,
                 const Poly_Impl::Sys<Gens>& gs);
bool
satisfies_all(const Topol topol, const Gen& g,
              const Poly_Impl::Sys<Cons>& cs);
bool
satisfies_all(const Topol topol, const Gen& g,
              const Poly_Impl::Sys<Cons>& cs1,
              const Poly_Impl::Sys<Cons>& cs2);
bool
satisfy_all(const Topol topol,
            const Poly_Impl::Sys<Gens>& gs,
            const Poly_Impl::Sys<Cons>& cs);
bool
satisfy_all(const Topol topol,
            const Poly_Impl::Sys<Gens>& gs1,
            const Poly_Impl::Sys<Gens>& gs2,
            const Poly_Impl::Sys<Cons>& cs1,
            const Poly_Impl::Sys<Cons>& cs2);
} // namespace sp

namespace detail {

template <typename Rows>
void
gauss(Rows& sg_rows) {
  const auto num_sg = num_rows(sg_rows);
  if (num_sg == 0)
    return;
  auto max_dim = sg_rows[0].linear_expr().last_nonzero();
  for (const auto& sg : sg_rows) {
    auto dim = sg.linear_expr().last_nonzero();
    if (dim > max_dim)
      max_dim = dim;
  }

  dim_type rank = 0;
  for (auto j : bwd_range(max_dim + 1)) {
    auto v_j = Var(j);
    for (auto i : range(rank, num_sg)) {
      auto& sg_i = sg_rows[i];
      assert(is_singular(sg_i));
      // Search for a non-zero coefficient in the j-th column.
      if (sg_i.coeff(v_j) == 0)
        continue;
      // Pivot found: make it the rank-th row.
      auto& sg_r = sg_rows[rank];
      if (i != rank) {
        using std::swap;
        swap(sg_i, sg_r);
      }
      // Eliminate v_j from other rows.
      for (auto k : range(i + 1, num_sg)) {
        auto& sg_k = sg_rows[k];
        if (sg_k.coeff(v_j) != 0)
          sg_k.linear_combine(sg_r, j);
      }
      ++rank;
      break;
    }
  }
  sg_rows.resize(rank);
}

template <typename Rows>
void
back_substitute(Rows& sg_rows, Rows& sk_rows) {
#ifndef NDEBUG
  using Row = typename Rows::value_type;
  auto is_sing = [](const Row& r) { return is_singular(r); };
  assert(all_of(sg_rows, is_sing));
  assert(none_of(sk_rows, is_sing));
#endif

  for (auto k : bwd_index_range(sg_rows)) {
    auto& sg_k = sg_rows[k];
    auto j = sg_k.linear_expr().last_nonzero();
    auto v_j = Var(j);
    // Go through sing rows above `sg_k'.
    for (auto i : bwd_range(k)) {
      auto& sg_i = sg_rows[i];
      if (sg_i.coeff(v_j) != 0)
        sg_i.linear_combine(sg_k, j);
    }
    // Force `sg_k[j]' to be positive.
    const bool have_to_negate = (sg_k.coeff(v_j) < 0);
    if (have_to_negate) {
      neg_assign(sg_k.impl().expr);
      neg_assign(sg_k.impl().inhomo);
    }
    // Go through all sk rows.
    for (auto& sk : sk_rows) {
      if (sk.coeff(v_j) != 0)
        sk.linear_combine(sg_k, j);
    }
    // Restore strong normalization.
    if (have_to_negate) {
      neg_assign(sg_k.impl().expr);
      neg_assign(sg_k.impl().inhomo);
    }
    assert(sg_k.check_inv());
  }
}

template <typename Sys>
inline void
normalize_for_hash(const Sys& sys) {
  // Only logically const.
  Sys& s = const_cast<Sys&>(sys);
  gauss(s.sing_rows);
  back_substitute(s.sing_rows, s.sk_rows);
}

using Con_Sys = Poly_Impl::Sys<Cons>;
using Gen_Sys = Poly_Impl::Sys<Gens>;

#if PPLITE_NOISY_CONVERSION

template <typename Dst_Sys>
void output_counters(const char* kind, dim_type k, const Dst_Sys& dst);

extern template void
output_counters<Con_Sys>(const char*, dim_type, const Con_Sys&);
extern template void
output_counters<Gen_Sys>(const char*, dim_type, const Gen_Sys&);

#define OUTPUT_COUNTERS_IF_NOISY(kind, k, dst) \
  { output_counters(kind, k, dst); }

#else // !PPLITE_NOISY_CONVERSION

#define OUTPUT_COUNTERS_IF_NOISY(kind, k, dst)

#endif // !PPLITE_NOISY_CONVERSION

// Only needed in gen-to-con minimize.
// Use maybe_ helpers to call from templatic contexts.
void ensure_strict_pos(Cons& sk_rows, NS_Rows& ns_rows);
inline void
maybe_ensure_strict_pos(Cons& sk_rows, NS_Rows& ns_rows) {
  ensure_strict_pos(sk_rows, ns_rows);
}
inline void
maybe_ensure_strict_pos(Gens&, NS_Rows&) {
  // Nothing to do.
}

/* Note: used for debugging only. */
bool
check_cutters(const NS_Rows& ns_dst, const Cons& sk_dst,
              const Gens& sk_src, const Sat& sat_g);
inline bool
maybe_check_cutters(const NS_Rows&, const Gens&, const Cons&, const Sat&) {
  // Nothing to check: only for gen-to-con conversion.
  return true;
}
inline bool
maybe_check_cutters(const NS_Rows& ns_dst, const Cons& sk_dst,
                    const Gens& sk_src, const Sat& sat_g) {
  return check_cutters(ns_dst, sk_dst, sk_src, sat_g);
}

/* Note: used for debugging only. */
bool
check_efc(const NS_Rows& cs_ns, const Cons& cs_sk,
          const Gens& gs_sk, const Sat& sat_g, const Sat& sat_c);
inline bool
maybe_check_efc(const NS_Rows&, const Gens&,
                const Cons&, const Sat&, const Sat&) {
  // Nothing to check: only for gen-to-con conversion.
  return true;
}
inline bool
maybe_check_efc(const NS_Rows& cs_ns, const Cons& cs_sk,
                const Gens& gs_sk, const Sat& sat_g, const Sat& sat_c) {
  return check_efc(cs_ns, cs_sk, gs_sk, sat_g, sat_c);
}

Index_Set
find_efc(const NS_Rows& cs_ns,
         const Cons& cs_sk, const Gens& gs_sk,
         const Sat& sat_g, const Sat& sat_c);

// Returns <idx1, idx2>, where:
//  * idx1 is the index of the positivity constraint in x.cs.sk_rows
//    or not_a_dim() if it is not there;
//  * idx2 is the index of the efc (empty face cutter) in x.cs.ns_rows
//    or not_a_dim() if it is not there.
inline std::pair<dim_type, dim_type>
get_pos_efc_info(const Poly_Impl& x) {
  assert(x.is_minimized() && !x.is_empty());
  Index_Set maybe_efc = find_efc(x.cs.ns_rows,
                                 x.cs.sk_rows, x.gs.sk_rows,
                                 x.sat_g, x.sat_c);
  auto sz = maybe_efc.size();
  if (sz == 1)
    // found the (skeleton) positivity constraint
    return { *maybe_efc.begin(), not_a_dim() };
  if (sz == 0)
    // efc is made redundant by other constraints
    return { not_a_dim(), not_a_dim() };

  // efc is (or should be) a non-skel constraint
  const auto& ns = x.cs.ns_rows;
  auto iter = std::find(ns.begin(), ns.end(), maybe_efc);
  if (iter == ns.end())
    return { not_a_dim(), not_a_dim() };
  else
    return { not_a_dim(), std::distance(ns.begin(), iter) };
}

/* Note: used for debugging only. */
bool check_redundant_ns_by_ns(const NS_Rows& ns_rows);

inline void
strict_ineqs_become_nonstrict_ineqs(Cons& cons) {
  for (auto& c : cons) {
    if (c.is_strict_inequality() && !c.is_tautological())
      c.set_type(Con::NONSTRICT_INEQUALITY);
  }
}
inline void
closure_points_become_points(Gens& gens) {
  for (auto& g : gens) {
    if (g.is_closure_point())
      g.set_type(Gen::POINT);
  }
}

template <typename Rows, typename Indices>
inline void
remove_rows(const Indices& indices,
            Rows& sk_rows, NS_Rows& ns_rows, Sat& sat) {
  erase_using_sorted_indices(sk_rows, indices);
  // Remap all supports to take into account the elements erased.
  for (auto& ns : ns_rows)
    ns.remove_all(indices);
  erase_using_sorted_indices(sat.impl().rows, indices);
}

inline void
make_singular(Con& c) {
  if (!c.is_equality()) {
    c.set_type(Con::EQUALITY);
    c.sign_normalize();
  }
}
inline void
make_singular(Gen& g) {
  if (!g.is_line()) {
    g.set_type(Gen::LINE);
    g.sign_normalize();
  }
}

template <typename SK_Rows>
inline dim_type
count_singular(const SK_Rows& sk_rows) {
  return std::count_if(sk_rows.begin(), sk_rows.end(), is_singular);
}

inline void
ensure_skel(Con& c) {
  if (c.is_strict_inequality())
    c.set_type(Con::NONSTRICT_INEQUALITY);
}
inline void
ensure_skel(Gen& g) {
  if (g.is_point())
    g.set_type(Gen::CLOSURE_POINT);
}

inline void
promote_singleton(Con& c) {
  assert(c.is_nonstrict_inequality());
  c.set_type(Con::STRICT_INEQUALITY);
}

inline void
promote_singleton(Gen& g) {
  assert(g.is_closure_point());
  g.set_type(Gen::POINT);
}

// For each ns in ns_rows, if it has only one supporter,
// change its supporter's type to point/strict-ineq in sk_rows
// and remove it from ns_rows.
template <typename Rows>
void
promote_singletons(NS_Rows& ns_rows, Rows& sk_rows) {
  Dims erasing;
  for (auto i : index_range(ns_rows)) {
    const auto& ns = ns_rows[i];
    if (ns.size() == 1) {
      dim_type supp_idx = *(ns.begin());
      promote_singleton(sk_rows[supp_idx]);
      erasing.push_back(i);
    }
  }
  erase_using_sorted_indices(ns_rows, erasing);
}

//////////////////////////////////////////////////////////
// Declarations of methods defined in Poly_conv.cc
// to be used in Poly_split.cc
//////////////////////////////////////////////////////////

// The possible results of a conversion iteration
// (note: EMPTY can only be returned by con-to-gen conversions).
enum class Conv_Iter_Result { OK, EMPTY, REDUNDANT };

void sk_partition(const Integers& sp, Ranges& sk_ranges);

// Returns true if we find in partition `part' of constraints `rows'
// at least two non-strict inequality constraints.
template <typename Part>
bool
can_have_nonskel(const Part& part, const Cons& rows) {
  bool seen_nonstrict = false;
  for (auto i : part) {
    if (rows[i].is_nonstrict_inequality()) {
      if (seen_nonstrict)
        return true;
      else
        seen_nonstrict = true;
    }
  }
  return false;
}

// Returns true if we find in partition `part' of generators `rows'
// both a closure point and another closure point or ray.
template <typename Part>
bool
can_have_nonskel(const Part& part,
                 const Gens& rows) {
  bool seen_cp = false;
  bool seen_ray = false;
  for (auto i : part) {
    if (rows[i].is_closure_point()) {
      if (seen_cp || seen_ray)
        return true;
      else
        seen_cp = true;
    } else if (!seen_ray && rows[i].is_ray()) {
      if (seen_cp)
        return true;
      else
        seen_ray = true;
    }
  }
  return false;
}

template <typename Index_Part>
inline void
update_last_sat_column(Sat& sat, const Index_Part& index_part) {
  assert(sat.num_cols() > 0);
  const auto last = sat.num_cols() - 1;
  for (auto i : index_part)
    sat[i].set(last);
}

// Classify non-skel row in NS=, NS+, NS-, NS+-
inline Range
where_is_ns(const Index_Set& ns,
            const Index_Set& sk_neg,
            const Index_Set& sk_pos) {
  assert(!ns.empty());
  const bool has_neg = ns.intersects(sk_neg);
  const bool has_pos = ns.intersects(sk_pos);
  return has_neg
    ? (has_pos ? Range::POS_NEG : Range::NEG)
    : (has_pos ? Range::POS : Range::EQ);
}

void sk_ranges_to_index_sets(const Ranges& sk_ranges,
                             Index_Set& iset_neg,
                             Index_Set& iset_eq,
                             Index_Set& iset_pos);
inline void
ns_partition(const NS_Rows& ns_rows,
             const Index_Set& iset_neg,
             const Index_Set& iset_pos,
             Ranges& ns_ranges) {
  ns_ranges.resize(ns_rows.size());
  for (auto i : bwd_index_range(ns_ranges))
    ns_ranges[i] = where_is_ns(ns_rows[i], iset_neg, iset_pos);
}

void sk_ranges_to_noneq_index_sets(const Ranges& sk_ranges,
                                   Index_Set& iset_neg,
                                   Index_Set& iset_pos);

Index_Set sk_ranges_to_eq_index_set(const Ranges& sk_ranges);

inline void
sk_to_ns_ranges(const Ranges& sk_ranges,
                const NS_Rows& ns_rows,
                Ranges& ns_ranges) {
  Index_Set iset_neg, iset_pos;
  sk_ranges_to_noneq_index_sets(sk_ranges, iset_neg, iset_pos);
  ns_partition(ns_rows, iset_neg, iset_pos, ns_ranges);
}

// When processing a skel-only element, project `supp' onto NS=;
// otherwise, when processing a strict_ineq or point, project it onto NS+.
inline void
supp_projection(Index_Set& supp, bool is_strict_inequality,
                const Index_Set& iset_neg, const Index_Set& iset_eq) {
  if (is_strict_inequality)
    supp -= iset_neg;
  else
    supp &= iset_eq;
}

// For each (skel) element indexed by `index_part',
// make sure it can not be a non-skeleton element; namely,
// map strict into non-strict ineqs and points into closure points.
template <typename Index_Part, typename SK_Rows>
inline void
points_become_closure_points(const Index_Part& index_part,
                             SK_Rows& sk_rows) {
  for (auto i : index_part)
    ensure_skel(sk_rows[i]);
}

template <typename Src_Row, typename Sing_Rows>
inline dim_type
check_sing_rows(const Src_Row& src, const Sing_Rows& rows) {
  using Sing_Row = typename Sing_Rows::value_type;
  auto iter = std::find_if(rows.cbegin(), rows.cend(),
                           [&src](const Sing_Row& sg) {
                             return sp::sign(sg, src) != 0;
                           });
  return std::distance(rows.cbegin(), iter);
}

inline Con::Type
combine_type(Con::Type c_pos, Con::Type c_neg,
             const bool is_point) {
  assert(c_pos != Con::EQUALITY);
  assert(c_neg != Con::EQUALITY);
  if (is_point)
    return Con::NONSTRICT_INEQUALITY;
  return (c_pos == Con::NONSTRICT_INEQUALITY)
    ? c_neg
    : Con::STRICT_INEQUALITY;
}

inline Gen::Type
combine_type(Gen::Type g_pos, Gen::Type g_neg,
             bool is_strict_inequality) {
  assert(g_pos != Gen::LINE);
  assert(g_neg != Gen::LINE);
  switch (g_pos) {
  case Gen::RAY:
    return (is_strict_inequality && g_neg == Gen::POINT)
      ? Gen::CLOSURE_POINT
      : g_neg;
  case Gen::POINT:
    return is_strict_inequality ? Gen::CLOSURE_POINT : Gen::POINT;
  case Gen::CLOSURE_POINT:
    return (!is_strict_inequality && g_neg == Gen::POINT)
      ? Gen::POINT
      : Gen::CLOSURE_POINT;
  case Gen::LINE:
  default:
    PPLITE_UNREACH;
    break;
  }
}

// Returns true if LeVerge's quick non-adjacency test succeded.
inline bool
quick_non_adj_test(const dim_type min_sat,
                   const dim_type max_sat,
                   const dim_type new_satrow_ones) {
  const dim_type num_common_satur = max_sat - new_satrow_ones;
  return (num_common_satur < min_sat);
}

// If either `sat[pos]' or `sat[neg]' has exactly one more zeroes
// than `new_satrow', then `dest_rows[pos]' and `dest_rows[neg]'
// are adjacent. Equivalently, adjacency holds if `new_satrow_ones'
// is equal to 1 plus the maximum of `|sat[pos]|' and `|sat[neg]|'.
// Note: the test is symmetric; names `pos' and `neg' are arbitrary.
inline bool
quick_adj_test(const dim_type pos,
               const dim_type neg,
               const dim_type new_satrow_ones,
               const Sat& sat) {
  const dim_type pos_num_ones = sat[pos].count_ones();
  const dim_type neg_num_ones = sat[neg].count_ones();
  const dim_type max_ones = std::max(pos_num_ones, neg_num_ones);
  return (max_ones + 1 == new_satrow_ones);
}

inline bool
quick_adj_test(const dim_type pos,
               const dim_type neg,
               const dim_type new_satrow_ones,
               const Dims& sat_ones) {
  const dim_type pos_num_ones = sat_ones[pos];
  const dim_type neg_num_ones = sat_ones[neg];
  const dim_type max_ones = std::max(pos_num_ones, neg_num_ones);
  return (max_ones + 1 == new_satrow_ones);
}

bool
combinatorial_adj_test(const dim_type first,
                       const dim_type last,
                       const dim_type pos,
                       const dim_type neg,
                       const Bits& new_satrow,
                       const Sat& sat_src);

template <typename Src_Row, typename Dst_Sys>
void
process_violating_singular(const dim_type ex, const Src_Row& src,
                           Dst_Sys& dst, Sat& sat_src,
                           Ranges& sk_ranges, Ranges& ns_ranges);

template <typename Src_Row, typename SK_Rows>
void
combine_with_ex_sing(const dim_type ex,
                     const Integer& sp_ex,
                     const Src_Row& src_k,
                     SK_Rows& sg_rows,
                     SK_Rows& sk_rows);

template <typename SK_Rows>
void
create_ns_with_ex_sing(SK_Rows& sk_rows,
                       NS_Rows& ns_rows,
                       const Ranges& sk_ranges,
                       const Ranges& ns_ranges,
                       dim_type ex);

// Return true if adjacency test succeded
inline bool
adj_test(const dim_type first,
         const dim_type last,
         const dim_type pos,
         const dim_type neg,
         const Bits& new_satrow,
         const dim_type new_satrow_ones,
         const dim_type min_sat,
         const dim_type max_sat,
         const Sat& sat_src) {
  // First, try the quick (but incomplete) non-adjacency test.
  if (quick_non_adj_test(min_sat, max_sat, new_satrow_ones))
    return false;
  // Then, try the quick (but incomplete) adjacency test.
  if (quick_adj_test(pos, neg, new_satrow_ones, sat_src))
    return true;
  // Finally, perform the complete (combinatorial) non-adjacency test.
  return combinatorial_adj_test(first, last, pos, neg, new_satrow, sat_src);
}

// Combines row_i with row_j editing row_i.
// Works for Row in { Con, Gen }.
template <typename Row>
inline void
combine_row(Row& row_i, const Row& row_j,
            const Integer& sp_i, const Integer& sp_j) {
  Linear_Expr::combine(row_i.impl().expr, row_i.impl().inhomo,
                       row_j.impl().expr, row_j.impl().inhomo,
                       sp_i, sp_j);
  row_i.strong_normalize();
}

inline bool
is_ex_singular(const Gen& g) { return g.is_ray(); }
inline bool
is_ex_singular(const Con& c) { return c.is_nonstrict_inequality(); }

inline bool
contains_a_point(const Gens& gs) { return has_point(gs); }
inline bool
contains_a_point(const Cons&) { return false; }

template <typename SK_Rows>
void
combine_sk_rows(Ranges& sk_ranges, SK_Rows& sk_rows,
                Integers& sp, Sat& sat_src, const bool is_strict_ineq,
                const dim_type min_sat, const dim_type max_sat);
void
add_to_minimal_set_of_ns(const Index_Set& new_ns, NS_Rows& set_of_ns);
void
add_to_minimal_set_of_ns(Index_Set&& new_ns, NS_Rows& set_of_ns);

void
remove_equivalent_generators(const Index_Set& F_g,
                             dim_type g,
                             const Bits& sat_F,
                             const Sat& sat_src,
                             Index_Set& search_range);

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
                               const bool is_strict_ineq);
template <typename SK_Rows, typename Index_Part>
inline void
create_ns_from_sk_with_face_enum(const SK_Rows& sk_rows,
                                 const Index_Part sk_part1,
                                 const Range range2,
                                 NS_Rows& set_of_ns,
                                 const Index_Set& iset_pos,
                                 const Index_Set& iset_eq,
                                 const Index_Set& iset_neg,
                                 const Sat& sat_src,
                                 const Sat& sat_dst,
                                 const bool is_strict_ineq) {
  for (auto i : sk_part1) {
    if (is_strict_ineq_or_point(sk_rows[i].type()))
      create_ns_with_enumerate_faces(Index_Set(i), range2,
                                     sk_rows, set_of_ns,
                                     iset_pos, iset_eq, iset_neg,
                                     sat_src, sat_dst, is_strict_ineq);
  }
}

template <typename SK_Rows, typename Index_Part>
inline void
create_ns_from_ns_with_face_enum(const NS_Rows& ns_rows,
                                 const SK_Rows& sk_rows,
                                 const Index_Part ns_part1,
                                 const Range range2,
                                 NS_Rows& set_of_ns,
                                 const Index_Set& iset_pos,
                                 const Index_Set& iset_eq,
                                 const Index_Set& iset_neg,
                                 const Sat& sat_src,
                                 const Sat& sat_dst,
                                 const bool is_strict_ineq) {
  for (auto i : ns_part1) {
    create_ns_with_enumerate_faces(ns_rows[i], range2,
                                   sk_rows, set_of_ns,
                                   iset_pos, iset_eq, iset_neg,
                                   sat_src, sat_dst, is_strict_ineq);
  }
}

// This works on both SK and NS for range1.
template <Range range1, typename SK_Rows>
inline void
create_ns_from_Q_with_face_enum(const Range range2,
                                const SK_Rows& sk_rows,
                                const NS_Rows& ns_rows,
                                const Ranges& sk_ranges,
                                const Ranges& ns_ranges,
                                NS_Rows& set_of_ns,
                                const Index_Set& iset_pos,
                                const Index_Set& iset_eq,
                                const Index_Set& iset_neg,
                                const Sat& sat_src,
                                const Sat& sat_dst,
                                const bool is_strict_ineq) {
  // Combine points in SK range1 with elements in SK range2.
  Index_Partition<range1> sk_part1(sk_ranges);
  create_ns_from_sk_with_face_enum(sk_rows, sk_part1, range2,
                                   set_of_ns,
                                   iset_pos, iset_eq, iset_neg,
                                   sat_src, sat_dst, is_strict_ineq);
  // Combine points in NS range1 with elements in SK range2.
  Index_Partition<range1> ns_part1(ns_ranges);
  create_ns_from_ns_with_face_enum(ns_rows, sk_rows, ns_part1, range2,
                                   set_of_ns,
                                   iset_pos, iset_eq, iset_neg,
                                   sat_src, sat_dst, is_strict_ineq);
}

template <typename SK_Rows>
void
create_ns_from_sk_with_adj_test(const SK_Rows& sk_rows,
                                const Ranges& sk_ranges,
                                NS_Rows& set_of_ns,
                                const Sat& sat_src,
                                const dim_type min_sat,
                                const dim_type max_sat);

template <typename SK_Rows>
inline void
create_ns_from_eq_points(const SK_Rows& sk_rows,
                         const NS_Rows& ns_rows,
                         const Ranges& sk_ranges,
                         const Ranges& ns_ranges,
                         NS_Rows& set_of_ns,
                         const Index_Set& iset_pos,
                         const Index_Set& iset_eq,
                         const Index_Set& iset_neg,
                         const Sat& sat_src,
                         const Sat& sat_dst,
                         const dim_type min_sat,
                         const dim_type max_sat) {
    // Combine points in SK= with elements in SK+.
    create_ns_from_sk_with_adj_test(sk_rows, sk_ranges, set_of_ns,
                                    sat_src, min_sat, max_sat);
    // Combine points in NS= with elements in SK+.
    Index_Partition<Range::EQ> ns_eq(ns_ranges);
    create_ns_from_ns_with_face_enum(ns_rows, sk_rows, ns_eq, Range::POS,
                                     set_of_ns,
                                     iset_pos, iset_eq, iset_neg,
                                     sat_src, sat_dst, true);
}

template<typename SK_Rows>
void
new_ns_from_set(const SK_Rows& sk_rows,
                NS_Rows& ns_dest,
                NS_Rows& set_of_ns,
                const Ranges& sk_ranges,
                const Ranges& ns_ranges,
                bool is_strict_ineq);

template <typename SK_Rows>
bool
check_redundant_ns_by_sk(const NS_Rows& ns_rows,
                         const SK_Rows& sk_rows);
bool
check_redundant_ns_by_ns(const NS_Rows& ns_rows);
bool
check_ns_rows(const NS_Rows& ns_rows);

inline void
make_non_singular(Gen& g, Integer& sp) {
  assert(g.is_line());
  g.set_type(Gen::RAY);
  if (sp < 0) {
    neg_assign(sp);
    neg_assign(g.linear_expr());
  }
  assert(g.check_inv());
}

inline void
make_non_singular(Con& c, Integer& sp) {
  assert(c.is_equality());
  c.set_type(Con::NONSTRICT_INEQUALITY);
  if (sp < 0) {
    neg_assign(sp);
    neg_assign(c.impl().expr);
    neg_assign(c.impl().inhomo);
  }
  assert(c.check_inv());
}

/***********************/

inline bool
has_no_points(const Cons&) {
  PPLITE_UNREACH;
  return true;
}
inline bool
has_no_points(const Gens& sk_rows) {
  return !has_point(sk_rows);
}

void
set_universe(dim_type space_dim, Con_Sys& cs, Gen_Sys& gs, Sat& sat_c);

inline void
init_dd(dim_type space_dim,
        Con_Sys& src, Con_Sys& /* pending */,
        Gen_Sys& dst, Sat& sat_c) {
  set_universe(space_dim, src, dst, sat_c);
}

void
init_dd(dim_type space_dim,
        Gen_Sys& src, Gen_Sys& pending,
        Con_Sys& dst, Sat& sat_g);

///////////////////////////////////////////////////////////////////////

// Function template: conversion.

// Template declaration.
template <bool con_to_gen, typename Src_Sys, typename Dst_Sys>
bool conversion(dim_type space_dim,
                Src_Sys& src, Src_Sys& pending,
                Dst_Sys& dst, Sat& sat_src);

// Declarations of explicit instantiations.
extern template
bool conversion<true, Con_Sys, Gen_Sys>(dim_type, Con_Sys&, Con_Sys&,
                                        Gen_Sys&, Sat&);
extern template
bool conversion<false, Gen_Sys, Con_Sys>(dim_type, Gen_Sys&, Gen_Sys&,
                                         Con_Sys&, Sat&);

///////////////////////////////////////////////////////////////////////

// Function template: simplify.

// Template declaration.
template <typename Src_Sys, typename Dst_Sys>
void simplify(dim_type space_dim, Src_Sys& src, Sat& sat_src, Sat& sat_dst,
              const Dst_Sys& dst);

// Declarations of explicit instantiations.
extern template
void
simplify<Con_Sys, Gen_Sys>(dim_type, Con_Sys&, Sat&, Sat&, const Gen_Sys&);
extern template
void
simplify<Gen_Sys, Con_Sys>(dim_type, Gen_Sys&, Sat&, Sat&, const Con_Sys&);

///////////////////////////////////////////////////////////////////////

template <bool con_to_gen, typename Src, typename Dst>
bool
add_and_minimize(const dim_type space_dim,
                 Src& src, Src& pending, Dst& dst,
                 Sat& sat_src, Sat& sat_dst) {
  const bool non_incr = src.empty();
#if PPLITE_CONVERSION_TIME_STATS
  LLOp_Clock clock(non_incr
                   ? (con_to_gen ? LLOp::CONV_C2G : LLOp::CONV_G2C)
                   : (con_to_gen ? LLOp::INCR_C2G : LLOp::INCR_G2C));
#endif // PPLITE_CONVERSION_TIME_STATS

  if (non_incr)
    init_dd(space_dim, src, pending, dst, sat_src);

  bool empty = conversion<con_to_gen>(space_dim, src, pending, dst, sat_src);

  if (empty)
    assert(con_to_gen);
  else
    simplify(space_dim, src, sat_src, sat_dst, dst);

  return empty;
}

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
             const Con& beta, Topol t);

/**
 * 2-way integral split method.
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
Z_split_poly(Poly_Impl& ph_then, Poly_Impl& ph_else,
             const Con& beta);

} // namespace detail

} // namespace pplite

#endif // !defined(pplite_Poly_min_hh)
