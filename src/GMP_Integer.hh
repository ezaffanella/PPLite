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

#ifndef pplite_GMP_Integer_hh
#define pplite_GMP_Integer_hh 1

#include "globals.hh"

#include <gmp.h>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

namespace pplite {

namespace detail {

inline size_t
hash(const mpz_t mp) {
  // A dumb, inefficient, just tentative hash function.
  constexpr unsigned long big_ui
    = 1UL << (std::numeric_limits<unsigned long>::digits - 1);
  return mpz_fdiv_ui(mp, big_ui);
}

} // namespace detail

class GMP_Integer {
private:
  mpz_t mp;
public:
  mpz_ptr impl() { return mp; }
  mpz_srcptr impl() const { return mp; }

  // Special members.
  GMP_Integer() noexcept { mpz_init(mp); }
  GMP_Integer(GMP_Integer const& x) { mpz_init_set(mp, x.mp); }
  void operator=(GMP_Integer const& x) { mpz_set(mp, x.mp); }
  GMP_Integer(GMP_Integer&& x) noexcept : GMP_Integer() { mpz_swap(mp, x.mp); }
  void operator=(GMP_Integer&& x) noexcept { mpz_swap(mp, x.mp); }
  ~GMP_Integer() { mpz_clear(mp); }

  // Non-special ctors (not explicit: behave as implicit conversions).
  GMP_Integer(signed int si) : GMP_Integer() { mpz_set_si(mp, si); }
  GMP_Integer(unsigned int ui) { mpz_init_set_ui(mp, ui); }
  GMP_Integer(signed long si) : GMP_Integer() { mpz_set_si(mp, si); }
  GMP_Integer(unsigned long ui) { mpz_init_set_ui(mp, ui); }

  // Explicit constructors from string.
  explicit GMP_Integer(const char* s) : GMP_Integer() {
    int res = mpz_set_str(mp, s, 10);
    assert(res == 0);
    if (res != 0) abort();
  }
  explicit GMP_Integer(const std::string& s) : GMP_Integer(s.c_str()) {}

  // Forbidden implicit conversions.
  GMP_Integer(bool) = delete;
  GMP_Integer(char) = delete;
  GMP_Integer(signed char) = delete;
  GMP_Integer(unsigned char) = delete;
  GMP_Integer(signed short) = delete;
  GMP_Integer(unsigned short) = delete;
  GMP_Integer(signed long long) = delete;
  GMP_Integer(unsigned long long) = delete;

  // Assignment operators.
  void operator=(signed int si) { mpz_set_si(mp, si); }
  void operator=(unsigned int ui) { mpz_set_ui(mp, ui); }
  void operator=(signed long si) { mpz_set_si(mp, si); }
  void operator=(unsigned long ui) { mpz_set_ui(mp, ui); }

  // Forbidden assignments.
  void operator=(bool) = delete;
  void operator=(char) = delete;
  void operator=(signed char) = delete;
  void operator=(unsigned char) = delete;
  void operator=(signed short) = delete;
  void operator=(unsigned short) = delete;
  void operator=(signed long long) = delete;
  void operator=(unsigned long long) = delete;

  // Conversion and assignment from GMP's mpz_t.
  GMP_Integer(mpz_t const z) { mpz_init_set(mp, z); }
  GMP_Integer& operator=(mpz_t const z) {
    mpz_set(mp, z);
    return *this;
  }

  void get_mpz(mpz_t dst) const { mpz_set(dst, mp); }

  size_t hash() const { return detail::hash(mp); }

  static const GMP_Integer& zero() {
    static PPLITE_TLS const GMP_Integer z;
    return z;
  }
  static const GMP_Integer& one() {
    static PPLITE_TLS const GMP_Integer z(1);
    return z;
  }

  bool is_zero() const { return mpz_sgn(mp) == 0; }
  bool is_one() const { return mpz_cmp_si(mp, 1) == 0; }
  bool is_pm1() const { return mpz_cmpabs_ui(mp, 1) == 0; }

  static void print(std::ostream& os, const mpz_t mp);
  static bool read(std::istream& os, mpz_t mp);

  void print(std::ostream& os) const { print(os, mp); }
  void print() const { print(std::cout, mp); }
  bool read(std::istream& is) { return read(is, mp); }

  void ascii_dump(std::ostream& s) const;
  bool ascii_load(std::istream& is);

}; // class GMP_Integer

inline size_t
external_memory_in_bytes(const __mpz_struct* ptr) {
  return ptr->_mp_alloc * PPLITE_SIZEOF_MP_LIMB_T;
}

inline bool operator==(GMP_Integer const& x, GMP_Integer const& y) {
  return 0 == mpz_cmp(x.impl(), y.impl());
}

inline bool operator==(GMP_Integer const& x, signed int si) {
  return 0 == mpz_cmp_si(x.impl(), si);
}

inline bool operator==(GMP_Integer const& x, signed long si) {
  return 0 == mpz_cmp_si(x.impl(), si);
}

inline bool operator==(signed int si, GMP_Integer const& x) {
  return 0 == mpz_cmp_si(x.impl(), si);
}

inline bool operator==(signed long si, GMP_Integer const& x) {
  return 0 == mpz_cmp_si(x.impl(), si);
}

inline bool operator==(GMP_Integer const& x, unsigned int ui) {
  return 0 == mpz_cmp_ui(x.impl(), ui);
}

inline bool operator==(GMP_Integer const& x, unsigned long ui) {
  return 0 == mpz_cmp_ui(x.impl(), ui);
}

inline bool operator==(unsigned int ui, GMP_Integer const& x) {
  return 0 == mpz_cmp_ui(x.impl(), ui);
}

inline bool operator==(unsigned long ui, GMP_Integer const& x) {
  return 0 == mpz_cmp_ui(x.impl(), ui);
}

inline bool operator!=(GMP_Integer const& x, GMP_Integer const& y) {
  return !(x == y);
}

inline bool operator!=(GMP_Integer const& x, signed int const& y) {
  return !(x == y);
}

inline bool operator!=(GMP_Integer const& x, signed long const& y) {
  return !(x == y);
}

inline bool operator!=(GMP_Integer const& x, unsigned int const& y) {
  return !(x == y);
}

inline bool operator!=(GMP_Integer const& x, unsigned long const& y) {
  return !(x == y);
}

inline bool operator<(GMP_Integer const& x, GMP_Integer const& y) {
  return (0 > mpz_cmp(x.impl(), y.impl()));
}

inline bool operator<(GMP_Integer const& x, signed int const& si) {
  return (0 > mpz_cmp_si(x.impl(), si));
}

inline bool operator<(GMP_Integer const& x, signed long const& si) {
  return (0 > mpz_cmp_si(x.impl(), si));
}

inline bool operator<(GMP_Integer const& x, unsigned int const& ui) {
  return (0 > mpz_cmp_ui(x.impl(), ui));
}

inline bool operator<(GMP_Integer const& x, unsigned long const& ui) {
  return (0 > mpz_cmp_ui(x.impl(), ui));
}

inline bool operator>(GMP_Integer const& x, GMP_Integer const& y) {
  return (y < x);
}

inline bool operator>(GMP_Integer const& x, signed int const& y) {
  return (y < x);
}

inline bool operator>(GMP_Integer const& x, signed long const& y) {
  return (y < x);
}

inline bool operator>(GMP_Integer const& x, unsigned int const& y) {
  return (y < x);
}

inline bool operator>(GMP_Integer const& x, unsigned long const& y) {
  return (y < x);
}

inline bool operator>=(GMP_Integer const& x, GMP_Integer const& y) {
  return !(x < y);
}

inline bool operator>=(GMP_Integer const& x, signed int const& y) {
  return !(x < y);
}

inline bool operator>=(GMP_Integer const& x, signed long const& y) {
  return !(x < y);
}

inline bool operator>=(GMP_Integer const& x, unsigned int const& y) {
  return !(x < y);
}

inline bool operator>=(GMP_Integer const& x, unsigned long const& y) {
  return !(x < y);
}

inline bool operator<=(GMP_Integer const& x, GMP_Integer const& y) {
  return !(x > y);
}

inline bool operator<=(GMP_Integer const& x, signed int const& y) {
  return !(x > y);
}

inline bool operator<=(GMP_Integer const& x, signed long const& y) {
  return !(x > y);
}

inline bool operator<=(GMP_Integer const& x, unsigned int const& y) {
  return !(x > y);
}

inline bool operator<=(GMP_Integer const& x, unsigned long const& y) {
  return !(x > y);
}

inline GMP_Integer& operator++(GMP_Integer& x) {
  mpz_add_ui(x.impl(), x.impl(), 1);
  return x;
}

inline GMP_Integer& operator--(GMP_Integer& x) {
  mpz_sub_ui(x.impl(), x.impl(), 1);
  return x;
}

inline void abs_assign(GMP_Integer& x) {
  mpz_abs(x.impl(), x.impl());
}

inline GMP_Integer abs(GMP_Integer const& x) {
  GMP_Integer res(x);
  abs_assign(res);
  return res;
}

inline void neg_assign(GMP_Integer& x) {
  mpz_neg(x.impl(), x.impl());
}

inline GMP_Integer neg(GMP_Integer const& x) {
  GMP_Integer res(x);
  neg_assign(res);
  return res;
}

inline GMP_Integer operator++(GMP_Integer& x, int) {
  GMP_Integer res(x);
  ++x;
  return res;
}

inline GMP_Integer operator--(GMP_Integer& x, int) {
  GMP_Integer res(x);
  --x;
  return res;
}

inline void swap(GMP_Integer& x, GMP_Integer& y) noexcept {
  mpz_swap(x.impl(), y.impl());
}

inline GMP_Integer&
operator+=(GMP_Integer& x, GMP_Integer const& y) {
  mpz_add(x.impl(), x.impl(), y.impl());
  return x;
}

inline GMP_Integer&
operator+=(GMP_Integer& x, signed int const& si) {
  x += GMP_Integer(si);
  return x;
}

inline GMP_Integer&
operator+=(GMP_Integer& x, signed long const& si) {
  x += GMP_Integer(si);
  return x;
}

inline GMP_Integer&
operator+=(GMP_Integer& x, unsigned int const& ui) {
  mpz_add_ui(x.impl(), x.impl(), ui);
  return x;
}

inline GMP_Integer&
operator+=(GMP_Integer& x, unsigned long const& ui) {
  mpz_add_ui(x.impl(), x.impl(), ui);
  return x;
}

inline GMP_Integer&
operator-=(GMP_Integer& x, GMP_Integer const& y) {
  mpz_sub(x.impl(), x.impl(), y.impl());
  return x;
}

inline GMP_Integer&
operator-=(GMP_Integer& x, signed int const& si) {
  x -= GMP_Integer(si);
  return x;
}

inline GMP_Integer&
operator-=(GMP_Integer& x, signed long const& si) {
  x -= GMP_Integer(si);
  return x;
}

inline GMP_Integer&
operator-=(GMP_Integer& x, unsigned int const& ui) {
  mpz_sub_ui(x.impl(), x.impl(), ui);
  return x;
}

inline GMP_Integer&
operator-=(GMP_Integer& x, unsigned long const& ui) {
  mpz_sub_ui(x.impl(), x.impl(), ui);
  return x;
}

inline GMP_Integer&
operator*=(GMP_Integer& x, GMP_Integer const& y) {
  mpz_mul(x.impl(), x.impl(), y.impl());
  return x;
}

inline GMP_Integer&
operator*=(GMP_Integer& x, signed int si) {
  mpz_mul_si(x.impl(), x.impl(), si);
  return x;
}

inline GMP_Integer&
operator*=(GMP_Integer& x, signed long si) {
  mpz_mul_si(x.impl(), x.impl(), si);
  return x;
}

inline GMP_Integer&
operator*=(GMP_Integer& x, unsigned int ui) {
  mpz_mul_ui(x.impl(), x.impl(), ui);
  return x;
}

inline GMP_Integer&
operator*=(GMP_Integer& x, unsigned long ui) {
  mpz_mul_ui(x.impl(), x.impl(), ui);
  return x;
}

inline GMP_Integer&
operator/=(GMP_Integer& x, GMP_Integer const& y) {
  assert(y != 0);
  mpz_tdiv_q(x.impl(), x.impl(), y.impl());
  return x;
}

inline GMP_Integer&
operator/=(GMP_Integer& x, signed int si) {
  x /= GMP_Integer(si);
  return x;
}

inline GMP_Integer&
operator/=(GMP_Integer& x, signed long si) {
  x /= GMP_Integer(si);
  return x;
}

inline GMP_Integer&
operator/=(GMP_Integer& x, unsigned int ui) {
  assert(ui != 0);
  mpz_tdiv_q_ui(x.impl(), x.impl(), ui);
  return x;
}

inline GMP_Integer&
operator/=(GMP_Integer& x, unsigned long ui) {
  assert(ui != 0);
  mpz_tdiv_q_ui(x.impl(), x.impl(), ui);
  return x;
}

inline GMP_Integer&
operator%=(GMP_Integer& x, GMP_Integer const& y) {
  assert(y != 0);
  mpz_tdiv_r(x.impl(), x.impl(), y.impl());
  return x;
}

inline GMP_Integer&
operator<<=(GMP_Integer& x, unsigned long ui) {
  mpz_mul_2exp(x.impl(), x.impl(), ui);
  return x;
}

inline GMP_Integer&
operator>>=(GMP_Integer& x, unsigned long ui) {
  mpz_tdiv_q_2exp(x.impl(), x.impl(), ui);
  return x;
}

inline GMP_Integer
operator+(GMP_Integer const& x, GMP_Integer const& y) {
  GMP_Integer res;
  mpz_add(res.impl(), x.impl(), y.impl());
  return res;
}

inline GMP_Integer
operator+(GMP_Integer const& x, signed int const& si) {
  return operator+(x, GMP_Integer(si));
}

inline GMP_Integer
operator+(GMP_Integer const& x, signed long const& si) {
  return operator+(x, GMP_Integer(si));
}

inline GMP_Integer
operator+(GMP_Integer const& x, unsigned int const& ui) {
  GMP_Integer res;
  mpz_add_ui(res.impl(), x.impl(), ui);
  return res;
}

inline GMP_Integer
operator+(GMP_Integer const& x, unsigned long const& ui) {
  GMP_Integer res;
  mpz_add_ui(res.impl(), x.impl(), ui);
  return res;
}

inline GMP_Integer
operator-(GMP_Integer const& x) {
  GMP_Integer res = x;
  neg_assign(res);
  return res;
}

inline GMP_Integer
operator-(GMP_Integer const& x, GMP_Integer const& y) {
  GMP_Integer res;
  mpz_sub(res.impl(), x.impl(), y.impl());
  return res;
}

inline GMP_Integer
operator-(GMP_Integer const& x, signed int const& si) {
  return operator-(x, GMP_Integer(si));
}

inline GMP_Integer
operator-(GMP_Integer const& x, signed long const& si) {
  return operator-(x, GMP_Integer(si));
}

inline GMP_Integer
operator-(GMP_Integer const& x, unsigned int const& ui) {
  GMP_Integer res;
  mpz_sub_ui(res.impl(), x.impl(), ui);
  return res;
}

inline GMP_Integer
operator-(GMP_Integer const& x, unsigned long const& ui) {
  GMP_Integer res;
  mpz_sub_ui(res.impl(), x.impl(), ui);
  return res;
}

inline GMP_Integer
operator*(GMP_Integer const& x, GMP_Integer const& y) {
  GMP_Integer res;
  mpz_mul(res.impl(), x.impl(), y.impl());
  return res;
}

inline GMP_Integer
operator*(GMP_Integer const& x, signed int si) {
  GMP_Integer res;
  mpz_mul_si(res.impl(), x.impl(), si);
  return res;
}

inline GMP_Integer
operator*(GMP_Integer const& x, signed long si) {
  GMP_Integer res;
  mpz_mul_si(res.impl(), x.impl(), si);
  return res;
}

inline GMP_Integer
operator*(GMP_Integer const& x, unsigned int ui) {
  GMP_Integer res;
  mpz_mul_ui(res.impl(), x.impl(), ui);
  return res;
}

inline GMP_Integer
operator*(GMP_Integer const& x, unsigned long ui) {
  GMP_Integer res;
  mpz_mul_ui(res.impl(), x.impl(), ui);
  return res;
}

inline GMP_Integer
operator/(GMP_Integer const& x, GMP_Integer const& y) {
  assert(y != 0);
  GMP_Integer res;
  mpz_tdiv_q(res.impl(), x.impl(), y.impl());
  return res;
}

inline GMP_Integer
operator/(GMP_Integer const& x, signed int si) {
  GMP_Integer res = x;
  res /= si;
  return res;
}

inline GMP_Integer
operator/(GMP_Integer const& x, signed long si) {
  GMP_Integer res = x;
  res /= si;
  return res;
}

inline GMP_Integer
operator/(GMP_Integer const& x, unsigned int ui) {
  GMP_Integer res = x;
  res /= ui;
  return res;
}

inline GMP_Integer
operator/(GMP_Integer const& x, unsigned long ui) {
  GMP_Integer res = x;
  res /= ui;
  return res;
}

inline GMP_Integer
operator%(GMP_Integer const& x, GMP_Integer const& y) {
  assert(y != 0);
  GMP_Integer res;
  mpz_tdiv_r(res.impl(), x.impl(), y.impl());
  return res;
}

inline GMP_Integer
operator<<(GMP_Integer const& x, unsigned long ui) {
  GMP_Integer res;
  mpz_mul_2exp(res.impl(), x.impl(), ui);
  return res;
}

inline GMP_Integer
operator>>(GMP_Integer const& x, unsigned long ui) {
  GMP_Integer res;
  mpz_tdiv_q_2exp(res.impl(), x.impl(), ui);
  return res;
}

inline int
compare(GMP_Integer const& x, GMP_Integer const& y) {
  return mpz_cmp(x.impl(), y.impl());
}

inline void
add_mul_assign(GMP_Integer& x,
               GMP_Integer const& y, GMP_Integer const& z) {
  mpz_addmul(x.impl(), y.impl(), z.impl());
}

inline void
sub_mul_assign(GMP_Integer& x,
               GMP_Integer const& y, GMP_Integer const& z) {
  mpz_submul(x.impl(), y.impl(), z.impl());
}

inline void
gcd_assign(GMP_Integer& x,
           GMP_Integer const& y, GMP_Integer const& z) {
  mpz_gcd(x.impl(), y.impl(), z.impl());
}

inline void
lcm_assign(GMP_Integer& x,
           GMP_Integer const& y, GMP_Integer const& z) {
  mpz_lcm(x.impl(), y.impl(), z.impl());
}

inline void
exact_div_assign(GMP_Integer& x,
                 GMP_Integer const& y, GMP_Integer const& z) {
  assert(z != 0);
  assert((y % z) == 0);
  mpz_divexact(x.impl(), y.impl(), z.impl());
}

inline int
sgn(GMP_Integer const& x) {
  return mpz_sgn(x.impl());
}

NOTHROW_DEFAULT_AND_MOVES(GMP_Integer);

} // namespace pplite

#endif // !defined(pplite_GMP_Integer_hh)
