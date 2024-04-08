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

#ifndef pplite_Con_hh
#define pplite_Con_hh 1

#include "globals.hh"
#include "Affine_Expr.hh"
#include "Integer.hh"
#include "Linear_Expr.hh"
#include "Var.hh"

#include <algorithm>
#include <cassert>
#include <functional>
#include <iostream>
#include <utility>
#include <vector>

namespace pplite {

class Con;

int compare(const Con& x, const Con& y);


class Con {
public:
  enum Type {
    EQUALITY,
    NONSTRICT_INEQUALITY,
    STRICT_INEQUALITY
  };

  struct Impl {
    Linear_Expr expr;
    Integer inhomo;
    Type type = EQUALITY;
    Impl() = default;
    Impl(Linear_Expr e, Integer i, Type t) noexcept
      : expr(std::move(e)), inhomo(std::move(i)), type(t) {}
  }; // struct Impl

private:
  Impl impl_;

public:
  Con() = default;
  Con(const Con& c) = default;
  Con(Con&& c) = default;
  Con& operator=(const Con& c) = default;
  Con& operator=(Con&& c) = default;
  ~Con() = default;

  Con(Linear_Expr expr, Integer inhomo, Type type)
    : impl_(std::move(expr), std::move(inhomo), type) {
    strong_normalize();
    assert(check_inv());
  }

  Con(Affine_Expr ae, Type type)
    : Con(std::move(ae.expr), std::move(ae.inhomo), type) {
  }

  dim_type space_dim() const { return linear_expr().space_dim(); }

  void set_space_dim(dim_type dim) { impl().expr.set_space_dim(dim); }

  void permute_space_dims_cycle(const Dims& cycle, dim_type d) {
    if (cycle.size() < 2)
      return;
    impl().expr.permute_space_dims_cycle(cycle, d);
    sign_normalize();
    assert(check_inv());
  }

  void shift_space_dims(dim_type start, dim_type n) {
    impl().expr.shift_space_dims(start, n);
  }

  Impl& impl() { return impl_; }
  const Impl& impl() const { return impl_; }

  Type type() const { return impl().type; };
  bool is_equality() const { return type() == EQUALITY; }
  bool is_inequality() const { return type() != EQUALITY; }
  bool is_nonstrict_inequality() const {
    return type() == NONSTRICT_INEQUALITY;
  }
  bool is_strict_inequality() const { return type() == STRICT_INEQUALITY; }

  Linear_Expr& linear_expr() { return impl().expr; }
  const Linear_Expr& linear_expr() const { return impl().expr; }
  const Integer& coeff(Var v) const { return impl().expr.get(v); }
  const Integer& inhomo_term() const { return impl().inhomo; }

  static Con zero_dim_false() {
    return Con(Linear_Expr(), Integer(1), EQUALITY);
  }
  static Con zero_dim_positivity() {
    return Con(Linear_Expr(), Integer(1), STRICT_INEQUALITY);
  }

  bool is_tautological() const;
  bool is_inconsistent() const;

  bool is_equal_to(const Con& y) const {
    return type() == y.type()
      && inhomo_term() == y.inhomo_term()
      && linear_expr().is_equal_to(y.linear_expr());
  }

  bool check_inv() const;

  void print(std::ostream& s) const;
  void print() const { print(std::cout); }

  void dump_type(std::ostream& s) const;
  void ascii_dump(std::ostream& s) const;

  bool load_type(std::istream& s);
  bool ascii_load(std::istream& s);

  void m_swap(Con& y) noexcept {
    using std::swap;
    swap(impl().expr, y.impl().expr);
    swap(impl().inhomo, y.impl().inhomo);
    swap(impl().type, y.impl().type);
  }

  void set_type(Type t) { impl().type = t; }

  bool is_line_or_equality() const { return is_equality(); }
  void set_is_line_or_equality() { set_type(EQUALITY); }

  void linear_combine(const Con& y, dim_type dim) {
    impl().expr.linear_combine(y.linear_expr(), dim,
                               impl().inhomo, y.inhomo_term());
    strong_normalize();
  }

  void sign_normalize() {
    if (is_equality())
      impl().expr.sign_normalize(impl().inhomo);
  }

  void strong_normalize() {
    impl().expr.normalize(impl().inhomo);
    sign_normalize();
  }

  bool check_strong_normalized() const {
    Con tmp = *this;
    tmp.strong_normalize();
    return compare(*this, tmp) == 0;
  }

}; // class Con

NOTHROW_DEFAULT_AND_MOVES(Con);

using Cons = std::vector<Con>;

namespace IO_Operators {

std::ostream& operator<<(std::ostream& s, const Con& c);
std::ostream& operator<<(std::ostream& s, const Con::Type& t);
std::ostream& operator<<(std::ostream& s, const Cons& cs);

} // namespace IO_Operators

inline bool
operator==(const Con& x, const Con& y) { return x.is_equal_to(y); }
inline bool
operator!=(const Con& x, const Con& y) { return !(x == y); }

inline int
compare(const Con& x, const Con& y) {
  const bool x_is_eq = x.is_equality();
  const bool y_is_eq = y.is_equality();
  if (x_is_eq != y_is_eq)
    // Equalities precede inequalities.
    return y_is_eq ? 1 : -1;
  const bool x_is_strict_ineq = x.is_strict_inequality();
  const bool y_is_strict_ineq = y.is_strict_inequality();
  if (x_is_strict_ineq != y_is_strict_ineq)
    // Strict inequalities follow non-strict inequalities.
    return x_is_strict_ineq ? 1 : -1;

  const int res = compare(x.linear_expr(), y.linear_expr());
  return (res == 0)
    ? compare(x.inhomo_term(), y.inhomo_term())
    : res;
}

inline Con
operator<(Linear_Expr e1, const Linear_Expr& e2) {
  neg_assign(e1);
  e1 += e2;
  return Con(std::move(e1), Integer::zero(), Con::STRICT_INEQUALITY);
}

inline Con
operator<(Var v1, Var v2) {
  return Con(v2 - v1, Integer::zero(), Con::STRICT_INEQUALITY);
}

inline Con
operator<(Linear_Expr e, Integer n) {
  neg_assign(e);
  return Con(std::move(e), std::move(n), Con::STRICT_INEQUALITY);
}

inline Con
operator<(Integer n, Linear_Expr e) {
  neg_assign(n);
  return Con(std::move(e), std::move(n), Con::STRICT_INEQUALITY);
}

inline Con
operator>(Linear_Expr e1, const Linear_Expr& e2) {
  e1 -= e2;
  return Con(std::move(e1), Integer::zero(), Con::STRICT_INEQUALITY);
}

inline Con
operator>(Var v1, Var v2) {
  return Con(v1 - v2, Integer::zero(), Con::STRICT_INEQUALITY);
}

inline Con
operator>(Linear_Expr e, Integer n) {
  neg_assign(n);
  return Con(std::move(e), std::move(n), Con::STRICT_INEQUALITY);
}

inline Con
operator>(Integer n, Linear_Expr e) {
  neg_assign(e);
  return Con(std::move(e), std::move(n), Con::STRICT_INEQUALITY);
}

inline Con
operator==(Linear_Expr e1, const Linear_Expr& e2) {
  e1 -= e2;
  return Con(std::move(e1), Integer::zero(), Con::EQUALITY);
}

inline Con
operator==(Var v1, Var v2) {
  return Con(v1 - v2, Integer::zero(), Con::EQUALITY);
}

inline Con
operator==(Linear_Expr e, Integer n) {
  neg_assign(n);
  return Con(std::move(e), std::move(n), Con::EQUALITY);
}

inline Con
operator==(Integer n, Linear_Expr e) {
  neg_assign(n);
  return Con(std::move(e), std::move(n), Con::EQUALITY);
}

inline Con
operator<=(Linear_Expr e1, const Linear_Expr& e2) {
  neg_assign(e1);
  e1 += e2;
  return Con(std::move(e1), Integer::zero(), Con::NONSTRICT_INEQUALITY);
}

inline Con
operator<=(Var v1, Var v2) {
  return Con(v2 - v1, Integer::zero(), Con::NONSTRICT_INEQUALITY);
}

inline Con
operator<=(Linear_Expr e, Integer n) {
  neg_assign(e);
  return Con(std::move(e), std::move(n), Con::NONSTRICT_INEQUALITY);
}

inline Con
operator<=(Integer n, Linear_Expr e) {
  neg_assign(n);
  return Con(std::move(e), std::move(n), Con::NONSTRICT_INEQUALITY);
}

inline Con
operator>=(Linear_Expr e1, const Linear_Expr& e2) {
  e1 -= e2;
  return Con(std::move(e1), Integer::zero(), Con::NONSTRICT_INEQUALITY);
}

inline Con
operator>=(Var v1, Var v2) {
  return Con(v1 - v2, Integer::zero(), Con::NONSTRICT_INEQUALITY);
}

inline Con
operator>=(Linear_Expr e, Integer n) {
  neg_assign(n);
  return Con(std::move(e), std::move(n), Con::NONSTRICT_INEQUALITY);
}

inline Con
operator>=(Integer n, Linear_Expr e) {
  neg_assign(e);
  return Con(std::move(e), std::move(n), Con::NONSTRICT_INEQUALITY);
}

inline Con
operator<(Affine_Expr e1, const Affine_Expr& e2) {
  neg_assign(e1);
  e1 += e2;
  return Con(std::move(e1), Con::STRICT_INEQUALITY);
}

inline Con
operator<(Affine_Expr e, const Integer& n) {
  neg_assign(e);
  e.inhomo += n;
  return Con(std::move(e), Con::STRICT_INEQUALITY);
}

inline Con
operator<(Affine_Expr e, Var v) {
  neg_assign(e);
  e.expr += v;
  return Con(std::move(e), Con::STRICT_INEQUALITY);
}

inline Con
operator<(const Integer& n, Affine_Expr e) {
  e.inhomo -= n;
  return Con(std::move(e), Con::STRICT_INEQUALITY);
}

inline Con
operator<(Var v, Affine_Expr e) {
  e.expr -= v;
  return Con(std::move(e), Con::STRICT_INEQUALITY);
}

inline Con
operator>(Affine_Expr e1, const Affine_Expr& e2) {
  e1 -= e2;
  return Con(std::move(e1), Con::STRICT_INEQUALITY);
}

inline Con
operator>(Affine_Expr e, const Integer& n) {
  e.inhomo -= n;
  return Con(std::move(e), Con::STRICT_INEQUALITY);
}

inline Con
operator>(Affine_Expr e, Var v) {
  e.expr -= v;
  return Con(std::move(e), Con::STRICT_INEQUALITY);
}

inline Con
operator>(const Integer& n, Affine_Expr e) {
  neg_assign(e);
  e.inhomo += n;
  return Con(std::move(e), Con::STRICT_INEQUALITY);
}

inline Con
operator>(Var v, Affine_Expr e) {
  neg_assign(e);
  e.expr += v;
  return Con(std::move(e), Con::STRICT_INEQUALITY);
}

inline Con
operator==(Affine_Expr e1, const Affine_Expr& e2) {
  e1 -= e2;
  return Con(std::move(e1), Con::EQUALITY);
}

inline Con
operator==(Affine_Expr e, const Integer& n) {
  e.inhomo -= n;
  return Con(std::move(e), Con::EQUALITY);
}

inline Con
operator==(Affine_Expr e, Var v) {
  e.expr -= v;
  return Con(std::move(e), Con::EQUALITY);
}

inline Con
operator==(const Integer& n, Affine_Expr e) {
  e.inhomo -= n;
  return Con(std::move(e), Con::EQUALITY);
}

inline Con
operator==(Var v, Affine_Expr e) {
  e.expr -= v;
  return Con(std::move(e), Con::EQUALITY);
}

inline Con
operator<=(Affine_Expr e1, const Affine_Expr& e2) {
  neg_assign(e1);
  e1 += e2;
  return Con(std::move(e1), Con::NONSTRICT_INEQUALITY);
}

inline Con
operator<=(Affine_Expr e, const Integer& n) {
  neg_assign(e);
  e.inhomo += n;
  return Con(std::move(e), Con::NONSTRICT_INEQUALITY);
}

inline Con
operator<=(Affine_Expr e, Var v) {
  neg_assign(e);
  e.expr += v;
  return Con(std::move(e), Con::NONSTRICT_INEQUALITY);
}

inline Con
operator<=(const Integer& n, Affine_Expr e) {
  e.inhomo -= n;
  return Con(std::move(e), Con::NONSTRICT_INEQUALITY);
}

inline Con
operator<=(Var v, Affine_Expr e) {
  e.expr -= v;
  return Con(std::move(e), Con::NONSTRICT_INEQUALITY);
}

inline Con
operator>=(Affine_Expr e1, const Affine_Expr& e2) {
  e1 -= e2;
  return Con(std::move(e1), Con::NONSTRICT_INEQUALITY);
}

inline Con
operator>=(Affine_Expr e, const Integer& n) {
  e.inhomo -= n;
  return Con(std::move(e), Con::NONSTRICT_INEQUALITY);
}

inline Con
operator>=(Affine_Expr e, Var v) {
  e.expr -= v;
  return Con(std::move(e), Con::NONSTRICT_INEQUALITY);
}

inline Con
operator>=(const Integer& n, Affine_Expr e) {
  neg_assign(e);
  e.inhomo += n;
  return Con(std::move(e), Con::NONSTRICT_INEQUALITY);
}

inline Con
operator>=(Var v, Affine_Expr e) {
  neg_assign(e);
  e.expr += v;
  return Con(std::move(e), Con::NONSTRICT_INEQUALITY);
}

inline void
swap(Con& x, Con& y) noexcept { x.m_swap(y); }

inline bool
has_strict_ineq(const Cons& cs) {
  return any_of(cs, std::mem_fn(&Con::is_strict_inequality));
}

inline bool
is_interval_con(const Con& c) {
  return is_interval_expr(c.linear_expr());
}

inline bool
is_proper_interval_con(const Con& c) {
  return is_proper_interval_expr(c.linear_expr());
}

namespace detail {

inline void
erase_higher_dims(Cons& cs, dim_type d) {
  for (auto& c : cs) {
    if (c.space_dim() > d) {
      c.set_space_dim(d);
      c.strong_normalize();
    }
  }
}

template <typename Indices>
Con materialize(const Indices& is, const Cons& cs) {
  Linear_Expr expr;
  Integer inhomo;
  for (auto i : is) {
    expr += cs[i].linear_expr();
    inhomo += cs[i].inhomo_term();
  }
  return Con(std::move(expr), std::move(inhomo), Con::STRICT_INEQUALITY);
}

// Note: this is kept in namespace detail, as it should only be called
// to remove space dims that are unconstrained, i.e., they always have
// a zero coefficient in all constraints.
template <typename Iter>
inline void
erase_space_dims(Cons& cs, Iter first, Iter last) {
  assert(first != last);
  for (auto& c : cs) {
    const auto dim = c.space_dim();
    Iter c_last = last;
    for ( ; first != c_last; --c_last) {
      Iter prev = c_last;
      --prev;
      if (*prev < dim)
        break;
    }
    erase_using_sorted_indices(c.linear_expr().impl(), first, c_last);
    assert(c.check_strong_normalized());
  }
}

// Build the "complement" of c.
// Note: it assumes c is not a trivial constraint. For equalities,
// the "complement" is only defined in the NNC case: it is the strict
// complement of the corresponding nonstrict inequality.
inline Con
complement_con(const Con& c, Topol t) {
  assert(!c.is_tautological() && !c.is_inconsistent());
  assert(t == Topol::NNC || c.is_nonstrict_inequality());
  auto type = (t == Topol::CLOSED || c.is_strict_inequality())
    ? Con::NONSTRICT_INEQUALITY
    : Con::STRICT_INEQUALITY;
  auto expr = c.linear_expr();
  auto inhomo = c.inhomo_term();
  neg_assign(inhomo);
  neg_assign(expr);
  return Con(std::move(expr), std::move(inhomo), type);
}

inline std::pair<Con, Con>
integral_complement_cons(const Con& c) {
  assert(not c.is_tautological() && not c.is_inconsistent());
  assert(not c.is_equality());
  auto expr_then = c.linear_expr();
  // Compute gcd and use to to refine inhomo_then.
  Integer expr_gcd = expr_then.gcd(0, c.space_dim());
  auto inhomo_then = c.inhomo_term();
  if (c.is_strict_inequality())
    --inhomo_then;
  inhomo_then -= abs(inhomo_then % expr_gcd);
  // Compute (rational non-strict) complement of refined constraint.
  auto expr_else = expr_then;
  neg_assign(expr_else);
  auto inhomo_else = inhomo_then;
  neg_assign(inhomo_else);
  // Refine inhomo_else.
  inhomo_else -= expr_gcd;
  Con c_then(std::move(expr_then), std::move(inhomo_then),
             Con::NONSTRICT_INEQUALITY);
  Con c_else(std::move(expr_else), std::move(inhomo_else),
             Con::NONSTRICT_INEQUALITY);
  return std::make_pair(std::move(c_then), std::move(c_else));
}

inline std::pair<Con, Con>
integral_complement_eq(const Con& c) {
  assert(not c.is_tautological() && not c.is_inconsistent());
  assert(c.is_equality());
  // Compute expr gcd (to refine inhomo terms).
  Integer expr_gcd = c.linear_expr().gcd(0, c.space_dim());
  // less than constraint needs to be negated.
  auto expr_lt = c.linear_expr();
  neg_assign(expr_lt);
  auto inhomo_lt = c.inhomo_term();
  neg_assign(inhomo_lt);
  // apply integral refinement
  --inhomo_lt;
  inhomo_lt -= abs(inhomo_lt % expr_gcd);
  Con c_lt(std::move(expr_lt), std::move(inhomo_lt), Con::NONSTRICT_INEQUALITY);
  // greater than constraint: no need to negate.
  const auto& expr_gt = c.linear_expr();
  auto inhomo_gt = c.inhomo_term();
  // apply integral refinement
  --inhomo_gt;
  inhomo_gt -= abs(inhomo_gt % expr_gcd);
  Con c_gt(expr_gt, std::move(inhomo_gt), Con::NONSTRICT_INEQUALITY);
  return std::make_pair(std::move(c_lt), std::move(c_gt));
}

inline bool
is_integral_inconsistent(const Con& c) {
  if (not c.is_equality())
    return c.is_inconsistent();
  // Compute gcd ignoring inhomo term.
  Integer expr_gcd = c.linear_expr().gcd(0, c.space_dim());
  auto modulo = c.inhomo_term() % expr_gcd;
  return not modulo.is_zero();
}

} // namespace detail

} // namespace pplite

#endif // !defined(pplite_Con_hh)
