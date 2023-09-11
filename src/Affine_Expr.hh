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

#ifndef pplite_Affine_Expr_hh
#define pplite_Affine_Expr_hh 1

#include "globals.hh"
#include "Integer.hh"
#include "Var.hh"
#include "Linear_Expr.hh"
#include <iostream>
#include <utility>

namespace pplite {

class Affine_Expr {
public: // Public access is meant.
  Integer inhomo;
  Linear_Expr expr;

  explicit Affine_Expr(Integer i = Integer()) noexcept
    : inhomo(std::move(i)), expr() {}
  explicit Affine_Expr(Linear_Expr e, Integer i = Integer())
    : inhomo(std::move(i)), expr(std::move(e)) {}

  void print(std::ostream& os) const;
  void print() const { print(std::cout); }

  Affine_Expr(Affine_Expr const& e) = default;
  Affine_Expr(Affine_Expr&& e) = default;
  Affine_Expr& operator=(Affine_Expr const& e) = default;
  Affine_Expr& operator=(Affine_Expr&& e) = default;
  ~Affine_Expr() = default;

  dim_type space_dim() const { return expr.space_dim(); }
  bool is_zero() const { return inhomo.is_zero() && expr.is_zero(); }

  void ascii_dump(std::ostream& s) const;
  bool ascii_load(std::istream& s);

  void m_swap(Affine_Expr& y) noexcept {
    using std::swap;
    swap(expr, y.expr);
    swap(inhomo, y.inhomo);
  }

  void normalize() { expr.normalize(inhomo); }
  void sign_normalize() { expr.sign_normalize(inhomo); }

}; // class Affine_Expr

inline Affine_Expr&
operator+=(Affine_Expr& a1, const Affine_Expr& a2) {
  a1.expr += a2.expr;
  a1.inhomo += a2.inhomo;
  return a1;
}
inline Affine_Expr&
operator+=(Affine_Expr& a1, const Linear_Expr& e2) {
  a1.expr += e2;
  return a1;
}
inline Affine_Expr&
operator+=(Affine_Expr& a, Var v) {
  a.expr += v;
  return a;
}
inline Affine_Expr&
operator+=(Affine_Expr& a, const Integer& n) {
  a.inhomo += n;
  return a;
}

inline Affine_Expr&
operator-=(Affine_Expr& a1, const Affine_Expr& a2) {
  a1.expr -= a2.expr;
  a1.inhomo -= a2.inhomo;
  return a1;
}
inline Affine_Expr&
operator-=(Affine_Expr& a1, const Linear_Expr& e2) {
  a1.expr -= e2;
  return a1;
}
inline Affine_Expr&
operator-=(Affine_Expr& a, Var v) {
  a.expr -= v;
  return a;
}
inline Affine_Expr&
operator-=(Affine_Expr& a, const Integer& n) {
  a.expr *= n;
  a.inhomo *= n;
  return a;
}

inline Affine_Expr&
operator*=(Affine_Expr& a, const Integer& n) {
  a.expr *= n;
  a.inhomo *= n;
  return a;
}

inline void
neg_assign(Affine_Expr& a) {
  neg_assign(a.expr);
  neg_assign(a.inhomo);
}

inline void
add_mul_assign(Affine_Expr& a, Integer const& n, Var v) {
  add_mul_assign(a.expr, n, v);
}
inline void
add_mul_assign(Affine_Expr& a, Integer const& n, const Linear_Expr& e) {
  add_mul_assign(a.expr, n, e);
}

inline void
sub_mul_assign(Affine_Expr& a, Integer const& n, Var v) {
  add_mul_assign(a.expr, n, v);
}
inline void
sub_mul_assign(Affine_Expr& a, Integer const& n, const Linear_Expr& e) {
  sub_mul_assign(a.expr, n, e);
}

inline Affine_Expr
operator+(Affine_Expr a1, const Affine_Expr& a2) { return a1 += a2; }
inline Affine_Expr
operator+(Affine_Expr a1, const Linear_Expr& e2) { return a1 += e2; }
inline Affine_Expr
operator+(Linear_Expr const& e1, Affine_Expr a2) { return a2 += e1; }
inline Affine_Expr
operator+(Affine_Expr a, const Integer& n) {
  a.inhomo += n;
  return a;
}
inline Affine_Expr
operator+(const Integer& n, Affine_Expr a) { return std::move(a) + n; }
inline Affine_Expr
operator+(Linear_Expr const& e, const Integer& n) { return Affine_Expr(e, n); }
inline Affine_Expr
operator+(const Integer& n, Linear_Expr const& e) { return Affine_Expr(e, n); }
inline Affine_Expr
operator+(Affine_Expr a, Var v) { return a += v; }
inline Affine_Expr
operator+(Var v, Affine_Expr a) { return a += v; }

// Binary minus.
inline Affine_Expr
operator-(Affine_Expr a1, const Affine_Expr& a2) { return a1 -= a2; }
inline Affine_Expr
operator-(Affine_Expr a1, const Linear_Expr& e2) { return a1 -= e2; }
inline Affine_Expr
operator-(Linear_Expr const& e1, Affine_Expr a2) {
  neg_assign(a2);
  return a2 += e1;
}
inline Affine_Expr
operator-(Affine_Expr a, const Integer& n) {
  a.inhomo -= n;
  return a;
}
inline Affine_Expr
operator-(const Integer& n, Affine_Expr a) {
  neg_assign(a);
  return a + n;
}
inline Affine_Expr
operator-(Linear_Expr const& e, const Integer& n) {
  return Affine_Expr(e, -n);
}
inline Affine_Expr
operator-(const Integer& n, Linear_Expr const& e) {
  return Affine_Expr(-e, n);
}
inline Affine_Expr
operator-(Affine_Expr a, Var v) { return a -= v; }
inline Affine_Expr
operator-(Var v, Affine_Expr a) {
  neg_assign(a);
  return a += v;
}

// Unary + and -.
inline Affine_Expr
operator+(Affine_Expr a) { return a; }
inline Affine_Expr
operator-(Affine_Expr a) {
  neg_assign(a);
  return a;
}

inline Affine_Expr
operator*(Integer const& n, Affine_Expr a) { return a *= n; }
inline Affine_Expr
operator*(Affine_Expr a, Integer const& n) { return n * std::move(a); }

namespace IO_Operators {

inline std::ostream&
operator<<(std::ostream& s, const Affine_Expr& a) {
  a.print(s);
  return s;
}

} // namespace IO_Operators

inline void swap(Affine_Expr& x, Affine_Expr& y) noexcept {
  x.m_swap(y);
}

NOTHROW_DEFAULT_AND_MOVES(Affine_Expr);

} // namespace pplite

#endif // !defined(pplite_Affine_Expr_hh)
