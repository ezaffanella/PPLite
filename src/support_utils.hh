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

#ifndef pplite_support_utils_hh
#define pplite_support_utils_hh 1

#include "utils.hh"
#include "Con.hh"
#include "Gen.hh"
#include "Bits.hh"
#include "Sat.hh"

namespace pplite {

namespace detail {

inline bool is_singular(Con::Type t) { return t == Con::EQUALITY; }
inline bool is_singular(Gen::Type t) { return t == Gen::LINE; }
inline bool is_singular(const Con& c) { return is_singular(c.type()); }
inline bool is_singular(const Gen& g) { return is_singular(g.type()); }

inline bool
is_strict_ineq_or_point(Con::Type t) { return t == Con::STRICT_INEQUALITY; }
inline bool
is_strict_ineq_or_point(Gen::Type t) { return t == Gen::POINT; }

inline bool
is_nonstrict_ineq_or_closure_point(Con::Type t) {
  return t == Con::NONSTRICT_INEQUALITY;
}
inline bool
is_nonstrict_ineq_or_closure_point(Gen::Type t) {
  return t == Gen::CLOSURE_POINT;
}

// Note: suffixes _g and _c are just for clarity.
// The function is meant to work also the other way round.
inline Bits
sat_all(const Index_Set& iset_g, const Sat& sat_c) {
  Bits row;
  for (auto i_g : iset_g)
    row |= sat_c[i_g];
  return row;
}

inline void
support_closure(Index_Set& iset, const Sat& sat_src, const Sat& sat_dst) {
  // Note: saturators are recorded as zeros.
  // Hence, in order to intersect them, we compute unions;
  // also, to map them to index sets, we use set complement.
  iset = sat_all(iset, sat_src);
  iset.complement_until(sat_src.num_cols());
  iset = sat_all(iset, sat_dst);
  iset.complement_until(sat_dst.num_cols());
}

// Returns true if it is maximal
// or if its saturated by rays only.
inline bool
is_empty_face_cutter(const Index_Set& supp,
                     const Sat& sat_g,
                     const Gens& gen_sk) {
  assert(sat_g.num_rows() > 0);
  if (supp.count_ones() == sat_g.num_rows())
    return true;

  Bits sat_row = sat_all(supp, sat_g);
  sat_row.complement_until(sat_g.num_cols());
  Index_Set saturators(std::move(sat_row));
  return all_of(saturators, [&gen_sk](dim_type i)
                            { return gen_sk[i].is_ray(); });
}

template <typename Rows>
inline bool
has_efc(const NS_Rows& ns_rows,
        const Sat& sat_g,
        const Rows& sk_rows) {
  return any_of(ns_rows, [&sat_g, &sk_rows](const Index_Set& ns) {
                          return is_empty_face_cutter(ns, sat_g, sk_rows);
                         });
}

} // namespace detail

} // namespace pplite

#endif // !defined(pplite_support_utils_hh)
