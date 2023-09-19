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

#ifndef pplite_Itv_hh
#define pplite_Itv_hh 1

#include "globals.hh"
#include "utils.hh"
#include "Con.hh"
#include "Integer.hh"
#include "Rational.hh"
#include "Var.hh"

#include <cassert>
#include <iostream>
#include <vector>

namespace pplite {

// Topologically closed 1-dim (rational) interval.
struct Itv {
  enum Kind { UNIVERSE, L_BOUNDED, U_BOUNDED, LU_BOUNDED, EMPTY };
  Kind kind = UNIVERSE;
  Rational lb;
  Rational ub;

  explicit Itv(Spec_Elem s = Spec_Elem::UNIVERSE) noexcept
    : kind(s == Spec_Elem::UNIVERSE ? UNIVERSE : EMPTY) {}

  static Itv zero() {
    Itv res;
    res.set_zero();
    return res;
  }

  bool check_inv() const {
    switch (kind) {
    case UNIVERSE:
    case EMPTY:
      return lb.is_zero() && ub.is_zero();
    case L_BOUNDED:
      return ub.is_zero();
    case U_BOUNDED:
      return lb.is_zero();
    case LU_BOUNDED:
      return lb <= ub;
    }
    return false;
  }

  static const Itv& empty_itv() {
    static PPLITE_TLS Itv res(Spec_Elem::EMPTY);
    return res;
  }

  bool is_empty() const { return kind == EMPTY; }
  bool is_universe() const { return kind == UNIVERSE; }
  bool has_lb() const { return kind == L_BOUNDED || kind == LU_BOUNDED; }
  bool has_ub() const { return kind == U_BOUNDED || kind == LU_BOUNDED; }
  bool inf_lb() const { return not has_lb(); }
  bool inf_ub() const { return not has_ub(); }
  bool is_bounded() const { return kind == LU_BOUNDED || kind == EMPTY; }
  bool is_singleton() const { return kind == LU_BOUNDED && lb == ub; }
  bool is_zero() const {
    return kind == LU_BOUNDED && lb.is_zero() && ub.is_zero();
  }

  bool is_disjoint_from(const Itv& y) const {
    const auto& x = *this;
    return x.is_empty() || y.is_empty()
      || (x.has_lb() && y.has_ub() && x.lb > y.ub)
      || (x.has_ub() && y.has_lb() && x.ub < y.lb);
  }
  bool intersects(const Itv& y) const {
    return !is_disjoint_from(y);
  }

  size_t hash() const {
    size_t res = kind;
    hash_combine(res, lb.hash());
    hash_combine(res, ub.hash());
    return res;
  }

  Rational length() const {
    assert(check_inv() && is_bounded());
    return ub - lb;
  }

  dim_type num_min_cons() const {
    switch (kind) {
    case UNIVERSE:
      return 0;
    case L_BOUNDED:
    case U_BOUNDED:
    case EMPTY:
      return 1;
    case LU_BOUNDED:
    default:
      return (lb == ub) ? 1 : 2;
    }
  }

  dim_type num_rays() const {
    switch (kind) {
    case UNIVERSE:
      return 2;
    case L_BOUNDED:
    case U_BOUNDED:
      return 1;
    case LU_BOUNDED:
    case EMPTY:
    default:
      return 0;
    }
  }

  bool contains(const Itv& y) const {
    const auto& x = *this;
    assert(x.check_inv());
    assert(y.check_inv());
    return (x.inf_lb() || (y.has_lb() && x.lb <= y.lb))
      && (x.inf_ub() || (y.has_ub() && x.ub >= y.ub));
  }

  bool contains(const Integer& num, const Integer& den) const {
    assert(check_inv());
    if (is_empty())
      return false;
    if (is_universe())
      return true;

    Rational r(num, den);
    if ((has_lb() && r < lb) || (has_ub() && r > ub))
      return false;
    return true;
  }

  bool operator==(const Itv& y) const {
    const auto& x = *this;
    assert(x.check_inv());
    assert(y.check_inv());
    return x.kind == y.kind
      && (x.inf_lb() || x.lb == y.lb)
      && (x.inf_ub() || x.ub == y.ub);
  }

  // Note: x < y implies !x.contains(y)
  bool operator<(const Itv& y) const {
    const auto& x = *this;
    if (y.inf_lb()) {
      if (x.has_lb())
        return true;
      if (y.inf_ub())
        return x.has_ub();
      if (x.inf_ub())
        return false;
      return x.ub < y.ub;
    }
    if (x.inf_lb())
      return false;
    assert(x.has_lb() && y.has_lb());
    switch (compare(x.lb, y.lb)) {
    case -1:
      return false;
    case 1:
      return true;
    case 0:
    default:
      if (x.inf_ub())
        return false;
      if (y.inf_ub())
        return true;
      assert(x.has_ub() && y.has_ub());
      return x.ub < y.ub;
    }
  }

  void set_empty() {
    kind = EMPTY;
    lb = Rational::zero();
    ub = Rational::zero();
  }
  void set_universe() {
    kind = UNIVERSE;
    lb = Rational::zero();
    ub = Rational::zero();
  }
  void set_zero() {
    kind = LU_BOUNDED;
    lb = Rational::zero();
    lb = Rational::zero();
  }

  void set_lb(Rational value) {
    kind = has_ub() ? LU_BOUNDED : L_BOUNDED;
    lb = std::move(value);
  }
  void set_ub(Rational value) {
    kind = has_lb() ? LU_BOUNDED : U_BOUNDED;
    ub = std::move(value);
  }
  void set_singleton(Rational value) {
    kind = LU_BOUNDED;
    lb = value;
    ub = std::move(value);
  }

  void unset_lb() {
    assert(!is_empty());
    if (has_lb()) {
      kind = is_bounded() ? U_BOUNDED : UNIVERSE;
      lb = Rational::zero();
    }
  }
  void unset_ub() {
    assert(!is_empty());
    if (has_ub()) {
      kind = is_bounded() ? L_BOUNDED : UNIVERSE;
      ub = Rational::zero();
    }
  }

  // Returns true if made empty
  bool glb_assign(const Itv& y) {
    auto& x = *this;
    if (y.has_lb()) {
      if (x.inf_lb())
        x.set_lb(y.lb);
      else if (x.lb < y.lb)
        x.lb = y.lb;
    }
    if (y.has_ub()) {
      if (x.inf_ub())
        x.set_ub(y.ub);
      else if (x.ub > y.ub)
        x.ub = y.ub;
    }
    // Check for inconsistency
    if (x.is_bounded() && (x.ub < x.lb)) {
      set_empty();
      return true;
    }
    return false;
  }

  void lub_assign(const Itv& y) {
    auto& x = *this;
    if (y.is_empty())
      return;
    if (x.is_empty()) {
      x = y;
      return;
    }
    if (x.has_lb()) {
      if (y.inf_lb())
        x.unset_lb();
      else if (x.lb > y.lb)
        x.lb = y.lb;
    }
    if (x.has_ub()) {
      if (y.inf_ub())
        x.unset_ub();
      else if (x.ub < y.ub)
        x.ub = y.ub;
    }
  }

  void widen_assign(const Itv& y) {
    auto& x = *this;
    // Adopt safe implementation (no inclusion precondition).
    if (y.is_empty() || (x == y))
      return;
    if (x.is_empty()) {
      x = y;
      return;
    }
    if (y.is_singleton()) {
      // Change of affine dimension.
      x.lub_assign(y);
      return;
    }
    if (x.has_lb() && y.has_lb() && x.lb < y.lb)
      x.unset_lb();
    if (x.has_ub() && y.has_ub() && x.ub > y.ub)
      x.unset_ub();
  }

  // Returns true if refinement makes it empty.
  bool refine_as_integral() {
    if (is_empty())
      return true;
    assert(not is_empty());
    if (has_lb()) lb.round_up();
    if (has_ub()) ub.round_down();
    if (is_bounded() && (ub < lb)) {
      set_empty();
      return true;
    }
    return false;
  }

  void complement_assign() {
    using std::swap;
    switch (kind) {
    case UNIVERSE:
      set_empty();
      return;
    case L_BOUNDED:
      kind = U_BOUNDED;
      swap(lb, ub);
      return;
    case U_BOUNDED:
      kind = L_BOUNDED;
      swap(lb, ub);
      return;
    case LU_BOUNDED:
      swap(lb, ub);
      return;
    case EMPTY:
      set_universe();
      return;
    }
  }

  void add_assign(const Itv& y) {
    if (is_empty())
      return;
    if (y.is_empty()) {
      set_empty();
      return;
    }
    if (is_universe())
      return;
    if (has_lb()) {
      if (y.has_lb())
        lb += y.lb;
      else
        unset_lb();
    }
    if (has_ub()) {
      if (y.has_ub())
        ub += y.ub;
      else
        unset_ub();
    }
  }

  void mul_assign(const Rational& r) {
    if (is_empty() || is_universe())
      return;
    if (has_lb())
      lb *= r;
    if (has_ub())
      ub *= r;
    if (sgn(r) < 0) {
      using std::swap;
      swap(lb, ub);
      if (kind == L_BOUNDED)
        kind = U_BOUNDED;
      else if (kind == U_BOUNDED)
        kind = L_BOUNDED;
    }
  }

  void print(std::ostream& os) const {
    using namespace IO_Operators;
    switch (kind) {
    case UNIVERSE:
      os << "(-inf, +inf)";
      break;
    case L_BOUNDED:
      os << "[" << lb << ", +inf)";
      break;
    case U_BOUNDED:
      os << "(-inf, " << ub << "]";
      break;
    case LU_BOUNDED:
      os << "[" << lb << ", " << ub << "]";
      break;
    case EMPTY:
      os << "{ }";
      break;
    }
  }

}; // struct Itv

NOTHROW_DEFAULT_AND_MOVES(Itv);

using Itvs = std::vector<Itv>;
using Itv_Expr = std::pair<Vars, Itvs>;

// Builds an itv from a the kind and inhomo of a constraint
// (i.e., ignoring the affine expression).
inline Itv
itv_from_con_inhomo(const Con& c) {
  if (c.linear_expr().is_zero()) {
    // Not a proper constraint.
    if (c.is_inconsistent())
      return Itv(Spec_Elem::EMPTY);
    else
      return Itv(Spec_Elem::UNIVERSE);
  }
  auto bound = Rational(-c.inhomo_term());
  Itv res;
  if (c.is_equality())
    res.set_singleton(bound);
  else
    res.set_lb(bound);
  return res;
}

inline Itv
itv_from_itv_con(const Con& c) {
  assert(is_interval_con(c));
  Itv res;
  auto d = c.linear_expr().first_nonzero();
  if (d == c.space_dim()) {
    // Not a proper constraint.
    if (c.is_inconsistent())
      res.set_empty();
    return res;
  }
  const auto& num = c.inhomo_term();
  const auto& den = c.linear_expr().get(d);
  auto b = Rational(num, den);
  neg_assign(b);
  if (c.is_equality())
    res.set_singleton(b);
  else if (den > 0)
    res.set_lb(b);
  else
    res.set_ub(b);
  return res;
}

inline Itv
split_itv(Itv& itv, const Con& c, bool integral) {
  assert(is_proper_interval_con(c));
  if (not integral) {
    // rational (closed) split
    assert(not c.is_equality());
    auto c_else = detail::complement_con(c, Topol::CLOSED);
    auto res = itv;
    itv.glb_assign(itv_from_itv_con(c));
    res.glb_assign(itv_from_itv_con(c_else));
    return res;
  }
  if (c.is_equality()) {
    // integral split on equality
#if 1 // Let header files be c++11 compliant
    auto cc = detail::integral_complement_eq(c);
    const auto& c_lt = cc.first;
    const auto& c_gt = cc.second;
#else // Structured bindings require c++17
    auto [c_lt, c_gt] = detail::integral_complement_eq(c);
#endif
    auto res_lt = itv;
    auto res_gt = itv;
    itv.glb_assign(itv_from_itv_con(c));
    itv.refine_as_integral();
    res_lt.glb_assign(itv_from_itv_con(c_lt));
    res_gt.glb_assign(itv_from_itv_con(c_gt));
    res_lt.lub_assign(res_gt);
    return res_lt;
  }
  // integral split on inequality
#if 1 // Let header files be c++11 compliant
  auto cc = detail::integral_complement_cons(c);
  const auto& c_then = cc.first;
  const auto& c_else = cc.second;
#else // Structured bindings require c++17
  auto [c_then, c_else] = detail::integral_complement_cons(c);
#endif
  auto res = itv;
  itv.glb_assign(itv_from_itv_con(c_then));
  res.glb_assign(itv_from_itv_con(c_else));
  return res;
}

inline Con
get_lb_con(Var var, const Itv& itv) {
  assert(itv.has_lb());
  return itv.lb.get_den() * var >= itv.lb.get_num();
}

inline Con
get_ub_con(Var var, const Itv& itv) {
  assert(itv.has_ub());
  return itv.ub.get_den() * var <= itv.ub.get_num();
}

inline Con
get_eq_con(Var var, const Itv& itv) {
  assert(itv.is_singleton());
  return itv.lb.get_den() * var == itv.lb.get_num();
}

namespace IO_Operators {

inline std::ostream&
operator<<(std::ostream& os, const Itv& itv) {
  itv.print(os);
  return os;
}

} // namespace IO_Operators

} // namespace pplite

#endif // !defined(pplite_Itv_hh)
