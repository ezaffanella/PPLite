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

#ifndef pplite_GMP_Rational_hh
#define pplite_GMP_Rational_hh 1

#include "globals.hh"
#include "utils.hh"
#include "GMP_Integer.hh"

#include <gmp.h>
#include <gmpxx.h>
#include <cassert>
#include <iostream>

namespace pplite {

class GMP_Rational {
private:
  mpq_t mp;
public:
  mpq_ptr impl() { return mp; }
  mpq_srcptr impl() const { return mp; }

  // Special members.
  GMP_Rational() noexcept { mpq_init(mp); }
  GMP_Rational(GMP_Rational const& x)
    : GMP_Rational() { mpq_set(mp, x.mp); }
  void operator=(GMP_Rational const& x) { mpq_set(mp, x.mp); }
  GMP_Rational(GMP_Rational&& x) noexcept
    : GMP_Rational() { mpq_swap(mp, x.mp); }
  void operator=(GMP_Rational&& x) noexcept { mpq_swap(mp, x.mp); }
  ~GMP_Rational() { mpq_clear(mp); }

  // Non-special ctors (explicit).
  explicit GMP_Rational(signed int n)
    : GMP_Rational() { mpq_set_si(mp, n, 1); }
  explicit GMP_Rational(signed long n)
    : GMP_Rational() { mpq_set_si(mp, n, 1); }
  explicit GMP_Rational(const GMP_Integer& n)
    : GMP_Rational() { mpq_set_z(mp, n.impl()); }
  GMP_Rational(const GMP_Integer& n, const GMP_Integer& d) {
    assert(d != 0);
    mpz_init_set(mpq_numref(mp), n.impl());
    mpz_init_set(mpq_denref(mp), d.impl());
    mpq_canonicalize(mp);
  }
  explicit GMP_Rational(double d)
    : GMP_Rational() { mpq_set_d(mp, d); }

  // Forbidden implicit conversions.
  GMP_Rational(bool) = delete;
  GMP_Rational(char) = delete;
  GMP_Rational(signed char) = delete;
  GMP_Rational(unsigned char) = delete;
  GMP_Rational(signed short) = delete;
  GMP_Rational(unsigned short) = delete;
  GMP_Rational(unsigned int) = delete;
  GMP_Rational(unsigned long) = delete;
  GMP_Rational(signed long long) = delete;
  GMP_Rational(unsigned long long) = delete;

  size_t hash() const {
    auto res = detail::hash(mpq_numref(mp));
    hash_combine(res, detail::hash(mpq_denref(mp)));
    return res;
  }

  static const GMP_Rational& zero() {
    static PPLITE_TLS const GMP_Rational q_zero;
    return q_zero;
  }
  static const GMP_Rational& one() {
    static PPLITE_TLS const GMP_Rational q_one(1);
    return q_one;
  }

  bool is_zero() const { return mpq_sgn(mp) == 0; }

  GMP_Integer get_num() const { return mpq_numref(mp); }
  GMP_Integer get_den() const { return mpq_denref(mp); }
  double get_double() const { return mpq_get_d(mp); }
  void get_mpq(mpq_t dst) const { mpq_set(dst, mp); }

  void round_up() {
    auto den = mpq_denref(mp);
    if (0 == mpz_cmp_si(den, 1))
      return;
    auto num = mpq_numref(mp);
    mpz_cdiv_q(num, num, den);
    mpz_set_si(den, 1);
  }

  void round_down() {
    auto den = mpq_denref(mp);
    if (0 == mpz_cmp_si(den, 1))
      return;
    auto num = mpq_numref(mp);
    mpz_fdiv_q(num, num, den);
    mpz_set_si(den, 1);
  }

  static void print(std::ostream& os, const mpq_t mp);
  static bool read(std::istream& os, mpq_t mp);

  void print(std::ostream& os) const { print(os, mp); }
  void print() const { print(std::cout, mp); }
  bool read(std::istream& is) { return read(is, mp); }

  void ascii_dump(std::ostream& s) const;
  bool ascii_load(std::istream& is);

}; // class GMP_Rational

inline size_t
external_memory_in_bytes(const __mpq_struct* ptr) {
  return external_memory_in_bytes(mpq_numref(ptr))
    + external_memory_in_bytes(mpq_denref(ptr));
}

inline bool operator==(GMP_Rational const& x, GMP_Rational const& y) {
  return mpq_equal(x.impl(), y.impl());
}

inline bool operator!=(GMP_Rational const& x, GMP_Rational const& y) {
  return !(x == y);
}

inline bool operator<(GMP_Rational const& x, GMP_Rational const& y) {
  return (0 > mpq_cmp(x.impl(), y.impl()));
}

inline bool operator>(GMP_Rational const& x, GMP_Rational const& y) {
  return (y < x);
}

inline bool operator>=(GMP_Rational const& x, GMP_Rational const& y) {
  return !(x < y);
}

inline bool operator<=(GMP_Rational const& x, GMP_Rational const& y) {
  return !(x > y);
}

inline void abs_assign(GMP_Rational& x) {
  mpq_abs(x.impl(), x.impl());
}

inline GMP_Rational abs(GMP_Rational const& x) {
  GMP_Rational res = x;
  abs_assign(res);
  return res;
}

inline void neg_assign(GMP_Rational& x) {
  mpq_neg(x.impl(), x.impl());
}

inline GMP_Rational neg(GMP_Rational const& x) {
  GMP_Rational res = x;
  neg_assign(res);
  return res;
}

inline void swap(GMP_Rational& x, GMP_Rational& y) noexcept {
  mpq_swap(x.impl(), y.impl());
}

inline GMP_Rational&
operator+=(GMP_Rational& x, GMP_Rational const& y) {
  mpq_add(x.impl(), x.impl(), y.impl());
  return x;
}

inline GMP_Rational&
operator-=(GMP_Rational& x, GMP_Rational const& y) {
  mpq_sub(x.impl(), x.impl(), y.impl());
  return x;
}

inline GMP_Rational&
operator*=(GMP_Rational& x, GMP_Rational const& y) {
  mpq_mul(x.impl(), x.impl(), y.impl());
  return x;
}

inline GMP_Rational&
operator/=(GMP_Rational& x, GMP_Rational const& y) {
  assert(!y.is_zero());
  mpq_div(x.impl(), x.impl(), y.impl());
  return x;
}

inline GMP_Rational
operator+(GMP_Rational const& x, GMP_Rational const& y) {
  GMP_Rational res;
  mpq_add(res.impl(), x.impl(), y.impl());
  return res;
}

inline GMP_Rational
operator-(GMP_Rational const& x) {
  GMP_Rational res = x;
  neg_assign(res);
  return res;
}

inline GMP_Rational
operator-(GMP_Rational const& x, GMP_Rational const& y) {
  GMP_Rational res;
  mpq_sub(res.impl(), x.impl(), y.impl());
  return res;
}

inline GMP_Rational
operator*(GMP_Rational const& x, GMP_Rational const& y) {
  GMP_Rational res;
  mpq_mul(res.impl(), x.impl(), y.impl());
  return res;
}

inline GMP_Rational
operator/(GMP_Rational const& x, GMP_Rational const& y) {
  assert(!y.is_zero());
  GMP_Rational res = x;
  mpq_div(res.impl(), x.impl(), y.impl());
  return res;
}

inline int
compare(GMP_Rational const& x, GMP_Rational const& y) {
  return mpq_cmp(x.impl(), y.impl());
}

inline void
add_mul_assign(GMP_Rational& x,
               GMP_Rational const& y, GMP_Rational const& z) {

  mpq_t mul;
  mpq_mul(mul, y.impl(), z.impl());
  mpq_add(x.impl(), x.impl(), mul);
}

inline void
sub_mul_assign(GMP_Rational& x,
               GMP_Rational const& y, GMP_Rational const& z) {
  mpq_t mul;
  mpq_mul(mul, y.impl(), z.impl());
  mpq_sub(x.impl(), x.impl(), mul);
}

inline int
sgn(GMP_Rational const& x) {
  return mpq_sgn(x.impl());
}

NOTHROW_DEFAULT_AND_MOVES(GMP_Rational);

inline void
exact_div_assign(GMP_Integer& x,
                 const GMP_Integer& y, const GMP_Rational& z) {
  mpz_mul(x.impl(), y.impl(), mpq_numref(z.impl()));
  mpz_divexact(x.impl(), x.impl(), mpq_denref(z.impl()));
}

inline GMP_Integer
lcm_dens(const std::vector<GMP_Rational>& rs) {
  GMP_Integer lcm = 1;
  for (const auto& r : rs)
    mpz_lcm(lcm.impl(), lcm.impl(), mpq_denref(r.impl()));
  return lcm;
}

inline GMP_Rational
pow_si(GMP_Rational const& x, signed long si) {
  GMP_Rational res;
  if (si == 0)
    return res;
  if (si > 0)
    res = x;
  else
    mpq_inv(res.impl(), x.impl());
  unsigned long ui = std::abs(si);
  auto mp = res.impl();
  mpz_pow_ui(mpq_numref(mp), mpq_numref(mp), ui);
  mpz_pow_ui(mpq_denref(mp), mpq_denref(mp), ui);
  return res;
}

} // namespace pplite

#endif // !defined(pplite_GMP_Rational_hh)
