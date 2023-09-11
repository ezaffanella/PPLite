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

#ifndef pplite_FLINT_Rational_hh
#define pplite_FLINT_Rational_hh 1

#include "globals.hh"
#include "utils.hh"
#include "FLINT_Integer.hh"

#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <gmp.h>
#include <gmpxx.h>

#include <cassert>
#include <cstddef>
#include <iostream>
#include <vector>

namespace pplite {

class FLINT_Rational {
private:
  fmpq_t mp;
public:
  using Impl = fmpq_t;
  Impl& impl() { return mp; }
  Impl const& impl() const { return mp; }

  // Special members.
  FLINT_Rational() noexcept { fmpq_init(mp); }
  FLINT_Rational(FLINT_Rational const& x)
    : FLINT_Rational() { fmpq_set(mp, x.mp); }
  void operator=(FLINT_Rational const& x) { fmpq_set(mp, x.mp); }
  FLINT_Rational(FLINT_Rational&& x) noexcept
    : FLINT_Rational() { fmpq_swap(mp, x.mp); }
  void operator=(FLINT_Rational&& x) noexcept { fmpq_swap(mp, x.mp); }
  ~FLINT_Rational() { fmpq_clear(mp); }

  // Non-special ctors (explicit).
  explicit FLINT_Rational(signed int n)
    : FLINT_Rational() { fmpq_set_si(mp, n, 1); }
  explicit FLINT_Rational(signed long n)
    : FLINT_Rational() { fmpq_set_si(mp, n, 1); }
  explicit FLINT_Rational(const FLINT_Integer& n, const FLINT_Integer& d = 1)
    : FLINT_Rational() { fmpq_set_fmpz_frac(mp, n.impl(), d.impl()); }
  explicit FLINT_Rational(double d)
    : FLINT_Rational() {
    mpq_class mpq(d);
    fmpq_set_mpq(mp, mpq.get_mpq_t());
  }

  // Forbidden implicit conversions.
  FLINT_Rational(bool) = delete;
  FLINT_Rational(char) = delete;
  FLINT_Rational(signed char) = delete;
  FLINT_Rational(unsigned char) = delete;
  FLINT_Rational(signed short) = delete;
  FLINT_Rational(unsigned short) = delete;
  FLINT_Rational(unsigned int) = delete;
  FLINT_Rational(unsigned long) = delete;
  FLINT_Rational(signed long long) = delete;
  FLINT_Rational(unsigned long long) = delete;

  operator mpq_class() const {
    mpq_class res;
    fmpq_get_mpq(res.get_mpq_t(), impl());
    return res;
  }

  size_t hash() const {
    auto res = detail::hash(fmpq_numref(mp));
    hash_combine(res, detail::hash(fmpq_denref(mp)));
    return res;
  }

  static const FLINT_Rational& zero() {
    static PPLITE_TLS const FLINT_Rational q_zero;
    return q_zero;
  }
  static const FLINT_Rational& one() {
    static PPLITE_TLS const FLINT_Rational q_one(1);
    return q_one;
  }

  bool is_zero() const { return fmpq_is_zero(mp); }

  FLINT_Integer get_num() const { return fmpq_numref(mp); }
  FLINT_Integer get_den() const { return fmpq_denref(mp); }

  void round_up() {
    auto den = fmpq_denref(mp);
    if (fmpz_is_one(den))
      return;
    auto num = fmpq_numref(mp);
    fmpz_cdiv_q(num, num, den);
    fmpz_set_si(den, 1);
  }

  void round_down() {
    auto den = fmpq_denref(mp);
    if (fmpz_is_one(den))
      return;
    auto num = fmpq_numref(mp);
    fmpz_fdiv_q(num, num, den);
    fmpz_set_si(den, 1);
  }

  void print(std::ostream& os) const {
    FLINT_Integer::print(os, fmpq_numref(mp));
    if (fmpz_is_one(fmpq_denref(mp)))
      return;
    os << "/";
    FLINT_Integer::print(os, fmpq_denref(mp));
  }

  void print() const { print(std::cout); }

  void ascii_dump(std::ostream& s) const;
  bool ascii_load(std::istream& is);

}; // class FLINT_Rational

inline size_t
external_memory_in_bytes(const fmpq* ptr) {
  return external_memory_in_bytes(fmpq_numref(ptr))
    + external_memory_in_bytes(fmpq_denref(ptr));
}

inline bool operator==(FLINT_Rational const& x, FLINT_Rational const& y) {
  return fmpq_equal(x.impl(), y.impl());
}

inline bool operator!=(FLINT_Rational const& x, FLINT_Rational const& y) {
  return !(x == y);
}

inline bool operator<(FLINT_Rational const& x, FLINT_Rational const& y) {
  return (0 > fmpq_cmp(x.impl(), y.impl()));
}

inline bool operator>(FLINT_Rational const& x, FLINT_Rational const& y) {
  return (y < x);
}

inline bool operator>=(FLINT_Rational const& x, FLINT_Rational const& y) {
  return !(x < y);
}

inline bool operator<=(FLINT_Rational const& x, FLINT_Rational const& y) {
  return !(x > y);
}

inline void abs_assign(FLINT_Rational& x) {
  fmpq_abs(x.impl(), x.impl());
}

inline FLINT_Rational abs(FLINT_Rational const& x) {
  FLINT_Rational res = x;
  abs_assign(res);
  return res;
}

inline void neg_assign(FLINT_Rational& x) {
  fmpq_neg(x.impl(), x.impl());
}

inline FLINT_Rational neg(FLINT_Rational const& x) {
  FLINT_Rational res = x;
  neg_assign(res);
  return res;
}

inline void swap(FLINT_Rational& x, FLINT_Rational& y) noexcept {
  fmpq_swap(x.impl(), y.impl());
}

inline FLINT_Rational&
operator+=(FLINT_Rational& x, FLINT_Rational const& y) {
  fmpq_add(x.impl(), x.impl(), y.impl());
  return x;
}

inline FLINT_Rational&
operator-=(FLINT_Rational& x, FLINT_Rational const& y) {
  fmpq_sub(x.impl(), x.impl(), y.impl());
  return x;
}

inline FLINT_Rational&
operator*=(FLINT_Rational& x, FLINT_Rational const& y) {
  fmpq_mul(x.impl(), x.impl(), y.impl());
  return x;
}

inline FLINT_Rational&
operator/=(FLINT_Rational& x, FLINT_Rational const& y) {
  assert(!y.is_zero());
  fmpq_div(x.impl(), x.impl(), y.impl());
  return x;
}

inline FLINT_Rational
operator+(FLINT_Rational const& x, FLINT_Rational const& y) {
  FLINT_Rational res;
  fmpq_add(res.impl(), x.impl(), y.impl());
  return res;
}

inline FLINT_Rational
operator-(FLINT_Rational const& x) {
  FLINT_Rational res = x;
  neg_assign(res);
  return res;
}

inline FLINT_Rational
operator-(FLINT_Rational const& x, FLINT_Rational const& y) {
  FLINT_Rational res;
  fmpq_sub(res.impl(), x.impl(), y.impl());
  return res;
}

inline FLINT_Rational
operator*(FLINT_Rational const& x, FLINT_Rational const& y) {
  FLINT_Rational res;
  fmpq_mul(res.impl(), x.impl(), y.impl());
  return res;
}

inline FLINT_Rational
operator/(FLINT_Rational const& x, FLINT_Rational const& y) {
  assert(!y.is_zero());
  FLINT_Rational res;
  fmpq_div(res.impl(), x.impl(), y.impl());
  return res;
}

inline int
compare(FLINT_Rational const& x, FLINT_Rational const& y) {
  return fmpq_cmp(x.impl(), y.impl());
}

inline void
add_mul_assign(FLINT_Rational& x,
               FLINT_Rational const& y, FLINT_Rational const& z) {
  fmpq_addmul(x.impl(), y.impl(), z.impl());
}

inline void
sub_mul_assign(FLINT_Rational& x,
               FLINT_Rational const& y, FLINT_Rational const& z) {
  fmpq_submul(x.impl(), y.impl(), z.impl());
}

inline int
sgn(FLINT_Rational const& x) {
  return fmpq_sgn(x.impl());
}

NOTHROW_DEFAULT_AND_MOVES(FLINT_Rational);

inline void
exact_div_assign(FLINT_Integer& x,
                 const FLINT_Integer& y, const FLINT_Rational& z) {
  fmpz_mul(x.impl(), y.impl(), fmpq_numref(z.impl()));
  fmpz_divexact(x.impl(), x.impl(), fmpq_denref(z.impl()));
}

inline FLINT_Integer
lcm_dens(const std::vector<FLINT_Rational>& rs) {
  FLINT_Integer lcm = 1;
  for (const auto& r : rs)
    fmpz_lcm(lcm.impl(), lcm.impl(), fmpq_denref(r.impl()));
  return lcm;
}

inline FLINT_Rational
pow_si(FLINT_Rational const& x, signed long si) {
  FLINT_Rational res;
  fmpq_pow_si(res.impl(), x.impl(), si);
  return res;
}

} // namespace pplite

#endif // !defined(pplite_FLINT_Rational_hh)
