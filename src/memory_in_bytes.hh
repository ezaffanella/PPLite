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

#ifndef pplite_memory_in_bytes_hh
#define pplite_memory_in_bytes_hh 1

#include "globals.hh"
#include "Integer.hh"
#include "Rational.hh"
#include "Itv.hh"
#include "Var.hh"
#include "Linear_Expr.hh"
#include "Con.hh"
#include "Gen.hh"
#include "Bits.hh"
#include "Sat.hh"
#include "Poly.hh"
#include "B_Poly.hh"
#include "F_Poly.hh"
#include "Two_Poly.hh"
#include "U_Poly.hh"
#include "PolySet.hh"

namespace pplite {

template <typename T>
inline size_t
total_memory_in_bytes(const T& t) {
  return sizeof(T) + external_memory_in_bytes(t);
}

inline size_t external_memory_in_bytes(dim_type) { return 0; }
inline size_t external_memory_in_bytes(Var) { return 0; }

inline size_t
external_memory_in_bytes(const Integer& i) {
  return external_memory_in_bytes(i.impl());
}

inline size_t
external_memory_in_bytes(const Rational& r) {
  return external_memory_in_bytes(r.impl());
}

inline size_t
external_memory_in_bytes(const Itv& itv) {
  return external_memory_in_bytes(itv.lb)
    + external_memory_in_bytes(itv.ub);
}

template <typename T>
size_t
external_memory_in_bytes(const std::vector<T>& vt) {
  size_t res = sizeof(T) * vt.capacity();
  for (const auto& t : vt)
    res += external_memory_in_bytes(t);
  return res;
}

template <typename T>
size_t
external_memory_in_bytes(const std::list<T>& lt) {
  size_t res = 0;
  for (const auto& t : lt)
    // approx: size of t + two pointers
    res += total_memory_in_bytes(t) + 2 * sizeof(void*);
  return res;
}

inline size_t
external_memory_in_bytes(const Linear_Expr& le) {
  return external_memory_in_bytes(le.impl());
}

inline size_t
external_memory_in_bytes(const Con& c) {
  return external_memory_in_bytes(c.impl().expr)
    + external_memory_in_bytes(c.impl().inhomo);
}

inline size_t
external_memory_in_bytes(const Gen& g) {
  return external_memory_in_bytes(g.impl().expr)
    + external_memory_in_bytes(g.impl().inhomo);
}

inline size_t
external_memory_in_bytes(const Bits& bits) {
  return bits.impl().capacity() * Bits::word_size;
}

inline size_t
external_memory_in_bytes(const Sat& sat) {
  return external_memory_in_bytes(sat.impl().rows);
}

template <typename T>
inline size_t
external_memory_in_bytes(const Poly_Impl::Sys<T>& sys) {
  return external_memory_in_bytes(sys.sing_rows)
    + external_memory_in_bytes(sys.sk_rows)
    + external_memory_in_bytes(sys.ns_rows);
}

inline size_t
external_memory_in_bytes(const Poly& ph) {
  return external_memory_in_bytes(ph.impl().cs)
    + external_memory_in_bytes(ph.impl().gs)
    + external_memory_in_bytes(ph.impl().sat_c)
    + external_memory_in_bytes(ph.impl().sat_g)
    + external_memory_in_bytes(ph.impl().cs_pending)
    + external_memory_in_bytes(ph.impl().gs_pending);
}

inline size_t
external_memory_in_bytes(const F_Poly& ph) {
  const auto& x = ph.impl();
  return external_memory_in_bytes(x.itvs)
    + external_memory_in_bytes(x.blocks)
    + external_memory_in_bytes(x.factors);
}

template <typename PH>
inline size_t
external_memory_in_bytes(const B_Wrap<PH>& ph) {
  auto res = external_memory_in_bytes(ph.impl_poly());
  if (ph.has_valid_bbox()) {
    const auto itvs = ph.impl_bbox().itvs;
    res += external_memory_in_bytes(itvs);
  }
  return res;
}

template <typename PH>
inline size_t
external_memory_in_bytes(const Stats<PH>& ph) {
  return external_memory_in_bytes(ph.base());
}

template <typename PH1, typename PH2>
inline size_t
external_memory_in_bytes(const Two_Poly<PH1, PH2>& ph) {
  return external_memory_in_bytes(ph.ph1)
    + external_memory_in_bytes(ph.ph2);
}

template <typename PH>
inline size_t
external_memory_in_bytes(const PolySet<PH>& ph) {
  return external_memory_in_bytes(ph.seq());
}

} // namespace pplite

#endif // !defined(memory_in_bytes_hh)
