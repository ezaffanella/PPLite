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

#ifndef pplite_Scalar_Prod_hh
#define pplite_Scalar_Prod_hh 1

#include "globals.hh"
#include "support_utils.hh"
#include "Low_Level_Stats.hh"
#include "Bits.hh"
#include "Integer.hh"
#include "Linear_Expr.hh"
#include "Con.hh"
#include "Gen.hh"

namespace pplite {
namespace sp {

inline void
add_sp_assign(Integer& z, const Linear_Expr& e1, const Linear_Expr& e2) {
#if PPLITE_LOW_LEVEL_COUNTERS
  LLOp_Stats::bump_scalar_prod();
#endif
  auto min_dim = std::min(e1.space_dim(), e2.space_dim());
  for (auto i : range(min_dim))
    add_mul_assign(z, e1[i], e2[i]);
}

inline void
add_sp_assign(Integer& z, const Index_Set& e1_nonzeroes,
              const Linear_Expr& e1, const Linear_Expr& e2) {
#if PPLITE_LOW_LEVEL_COUNTERS
  LLOp_Stats::bump_scalar_prod();
#endif
  auto e2_dim = e2.space_dim();
  for (auto i : e1_nonzeroes) {
    assert(not e1[i].is_zero());
    if (i >= e2_dim) return;
    add_mul_assign(z, e1[i], e2[i]);
  }
}

inline void
assign(Integer& z,
       const Linear_Expr& e1, const Integer& z1,
       const Linear_Expr& e2, const Integer& z2) {
  z = z1 * z2;
  add_sp_assign(z, e1, e2);
}

inline void
assign(Integer& z, const Index_Set& e1_nonzeroes,
       const Linear_Expr& e1, const Integer& z1,
       const Linear_Expr& e2, const Integer& z2) {
  z = z1 * z2;
  add_sp_assign(z, e1_nonzeroes, e1, e2);
}

template <typename Row1, typename Row2>
inline void
assign(Integer& z, const Row1& row1, const Row2& row2) {
  assign(z,
         row1.impl().expr, row1.impl().inhomo,
         row2.impl().expr, row2.impl().inhomo);
}

template <typename Row1, typename Row2>
inline void
assign(Integer& z, const Index_Set& row1_nonzeroes,
       const Row1& row1, const Row2& row2) {
  assign(z, row1_nonzeroes,
         row1.impl().expr, row1.impl().inhomo,
         row2.impl().expr, row2.impl().inhomo);
}

template <typename Row>
inline bool
is_sparse(const Row& row, Index_Set& nonzeroes) {
  const auto& expr = row.linear_expr();
  const auto dim = expr.space_dim();
  const auto sp_factor = 3;
  auto nz_limit = 1 + dim / sp_factor;
  for (auto i : bwd_range(dim)) {
    if (expr[i].is_zero())
      continue;
    if (nz_limit == 0)
      return false;
    --nz_limit;
    nonzeroes.set(i);
  }
  return true;
}

template <typename Src_Row, typename Dst_Rows>
inline void
compute_prods(const Src_Row& src_row,
              const Dst_Rows& dst_rows,
              Integers& sp) {
  const auto dst_size = dst_rows.size();
  sp.resize(dst_size);
  Index_Set src_nonzeroes;
  if (is_sparse(src_row, src_nonzeroes)) {
    for (auto i : bwd_index_range(dst_rows))
      assign(sp[i], src_nonzeroes, src_row, dst_rows[i]);
  } else {
    for (auto i : bwd_index_range(dst_rows))
      assign(sp[i], src_row, dst_rows[i]);
  }
}

inline int
sign(const Index_Set& e1_nonzeroes,
     const Linear_Expr& e1, const Linear_Expr& e2) {
  Integer z;
  add_sp_assign(z, e1_nonzeroes, e1, e2);
  return sgn(z);
}

template <typename Row1, typename Row2>
inline int
sign(const Row1& row1, const Row2& row2) {
  static PPLITE_TLS Integer z;
  assign(z, row1, row2);
  return sgn(z);
}

template <typename Row1, typename Row2>
inline int
sign(const Index_Set& row1_nonzeroes,
     const Row1& row1, const Row2& row2) {
  static PPLITE_TLS Integer z;
  assign(z, row1_nonzeroes, row1, row2);
  return sgn(z);
}

inline bool
sign_violation(int s, const Con& c, const Gen& g) {
  if (g.is_line())
    return s != 0;
  switch (c.type()) {
  case Con::EQUALITY:
    return s != 0;
  case Con::NONSTRICT_INEQUALITY:
    return s < 0;
  case Con::STRICT_INEQUALITY:
  default:
    return g.is_point() ? (s <= 0) : (s < 0);
  }
}

inline bool
violated_by(const Con& c, const Gen& g) {
  return sign_violation(sign(c, g), c, g);
}

inline bool
violated_by(const Index_Set& c_nonzeroes, const Con& c, const Gen& g) {
  return sign_violation(sign(c_nonzeroes, c, g), c, g);
}

template <typename Src_Row, typename Dst_Rows>
bool
compute_sat_row(const Src_Row& src, const Dst_Rows& dst_rows,
                Bits& sat_row) {
  Index_Set src_nz;
  bool sparse = is_sparse(src, src_nz);
  for (auto i : bwd_index_range(dst_rows)) {
    int s = sparse
      ? sign(src_nz, src, dst_rows[i])
      : sign(src, dst_rows[i]);
    if (s > 0)
      sat_row.set(i);
    else if (s < 0)
      return false;
    else if (detail::is_strict_ineq_or_point(src.type())
             && detail::is_strict_ineq_or_point(dst_rows[i].type()))
      return false;
  }
  return true;
}

template <typename GS>
inline bool
satisfied_by_all(const Con& c, const GS& gs) {
  Index_Set c_nz;
  if (is_sparse(c, c_nz))
    return none_of(gs, [&](const Gen& g) {
                         return violated_by(c_nz, c, g);
                       });
  else
    return none_of(gs, [&](const Gen& g) {
                         return violated_by(c, g);
                       });
}

template <typename CS, typename GS>
inline bool
satisfied_by_all(const CS& cs, const GS& gs) {
  return all_of(cs, [&gs](const Con& c) {
                      return satisfied_by_all(c, gs);
                    });
}

template <typename CS>
inline bool
satisfies_all(const Gen& g, const CS& cs) {
  return none_of(cs, [&g](const Con& c) { return violated_by(c, g); });
}

template <typename GS, typename CS>
inline bool
satisfy_all(const GS& gs, const CS& cs, Sat& sat) {
  for (auto i : bwd_index_range(gs)) {
    if (!compute_sat_row(gs[i], cs, sat[i]))
      return false;
  }
  return true;
}

} // namespace sp
} // namespace pplite

#endif // !defined(pplite_Scalar_Prod_hh)
