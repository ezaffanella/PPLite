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

#ifndef pplite_FLINT_Integer_hh
#define pplite_FLINT_Integer_hh 1

#include "globals.hh"

#include <flint/fmpz.h>
#include <gmp.h>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

namespace pplite {

namespace detail {

inline size_t
hash(const fmpz_t mp) {
  // A dumb, inefficient, just tentative hash function.
  constexpr unsigned long big_ui
    = 1UL << (std::numeric_limits<unsigned long>::digits - 1);
  return fmpz_fdiv_ui(mp, big_ui);
}

} // namespace detail

class FLINT_Integer {
private:
  fmpz_t mp;
public:
  using fmpz_ptr = fmpz*;
  using fmpz_srcptr = const fmpz*;
  fmpz_ptr impl() { return mp; }
  fmpz_srcptr impl() const { return mp; }

  // Special members.
  FLINT_Integer() noexcept { fmpz_init(mp); }
  FLINT_Integer(FLINT_Integer const& x) { fmpz_init_set(mp, x.mp); }
  void operator=(FLINT_Integer const& x) { fmpz_set(mp, x.mp); }
  FLINT_Integer(FLINT_Integer&& x) noexcept
    : FLINT_Integer() { fmpz_swap(mp, x.mp); }
  void operator=(FLINT_Integer&& x) noexcept { fmpz_swap(mp, x.mp); }
  ~FLINT_Integer() { fmpz_clear(mp); }

  // Non-special ctors (not explicit: behave as implicit conversions).
  FLINT_Integer(signed int si) { fmpz_init_set_si(mp, si); }
  FLINT_Integer(unsigned int ui) { fmpz_init_set_ui(mp, ui); }
  FLINT_Integer(signed long si) { fmpz_init_set_si(mp, si); }
  FLINT_Integer(unsigned long ui) { fmpz_init_set_ui(mp, ui); }

  // Explicit constructors from string.
  explicit FLINT_Integer(const char* s) : FLINT_Integer() {
    int res = fmpz_set_str(mp, s, 10);
    assert(res == 0);
    if (res != 0) abort();
  }
  explicit FLINT_Integer(const std::string& s) : FLINT_Integer(s.c_str()) {}

  // Forbidden implicit conversions.
  FLINT_Integer(bool) = delete;
  FLINT_Integer(char) = delete;
  FLINT_Integer(signed char) = delete;
  FLINT_Integer(unsigned char) = delete;
  FLINT_Integer(signed short) = delete;
  FLINT_Integer(unsigned short) = delete;
  FLINT_Integer(signed long long si) = delete;
  FLINT_Integer(unsigned long long ui) = delete;

  // Assignment operators.
  void operator=(signed int si) { fmpz_set_si(mp, si); }
  void operator=(unsigned int ui) { fmpz_set_ui(mp, ui); }
  void operator=(signed long si) { fmpz_set_si(mp, si); }
  void operator=(unsigned long ui) { fmpz_set_ui(mp, ui); }

  // Forbidden assignments.
  void operator=(bool) = delete;
  void operator=(char) = delete;
  void operator=(signed char) = delete;
  void operator=(unsigned char) = delete;
  void operator=(signed short) = delete;
  void operator=(unsigned short) = delete;
  void operator=(signed long long) = delete;
  void operator=(unsigned long long) = delete;

  // Conversion from fmpz_t.
  FLINT_Integer(const fmpz_t z) { fmpz_init_set(mp, z); }

  // Conversion and assignment from GMP's mpz_t.
  FLINT_Integer(const mpz_t z)
    : FLINT_Integer() { fmpz_set_mpz(mp, z); }
  FLINT_Integer& operator=(const mpz_t z) {
    fmpz_set_mpz(mp, z);
    return *this;
  }

  size_t hash() const { return detail::hash(mp); }

  static const FLINT_Integer& zero() {
    static PPLITE_TLS const FLINT_Integer z;
    return z;
  }
  static const FLINT_Integer& one() {
    static PPLITE_TLS const FLINT_Integer z(1);
    return z;
  }

  bool is_zero() const { return fmpz_is_zero(mp); }
  bool is_one() const { return fmpz_is_one(mp); }
  bool is_pm1() const { return fmpz_is_pm1(mp); }

  static void print(std::ostream& os, const fmpz_t mp);
  static bool read(std::istream& os, fmpz_t mp);

  void print(std::ostream& os) const { print(os, mp); }
  void print() const { print(std::cout, mp); }
  bool read(std::istream& is) { return read(is, mp); }

  void ascii_dump(std::ostream& s) const;
  bool ascii_load(std::istream& is);

}; // class FLINT_Integer

inline size_t
external_memory_in_bytes(const fmpz* ptr) {
  const auto& d = *ptr;
  return COEFF_IS_MPZ(d)
    ? sizeof(mpz_t) + COEFF_TO_PTR(d)->_mp_alloc * PPLITE_SIZEOF_MP_LIMB_T
    : 0;
}

inline bool operator==(FLINT_Integer const& x, FLINT_Integer const& y) {
  return fmpz_equal(x.impl(), y.impl());
}

inline bool operator==(FLINT_Integer const& x, signed int si) {
  return fmpz_equal_si(x.impl(), si);
}

inline bool operator==(FLINT_Integer const& x, signed long si) {
  return fmpz_equal_si(x.impl(), si);
}

inline bool operator==(signed int si, FLINT_Integer const& x) {
  return fmpz_equal_si(x.impl(), si);
}

inline bool operator==(signed long si, FLINT_Integer const& x) {
  return fmpz_equal_si(x.impl(), si);
}

inline bool operator==(FLINT_Integer const& x, unsigned int ui) {
  return fmpz_equal_ui(x.impl(), ui);
}

inline bool operator==(FLINT_Integer const& x, unsigned long ui) {
  return fmpz_equal_ui(x.impl(), ui);
}

inline bool operator==(unsigned int ui, FLINT_Integer const& x) {
  return fmpz_equal_ui(x.impl(), ui);
}

inline bool operator==(unsigned long ui, FLINT_Integer const& x) {
  return fmpz_equal_ui(x.impl(), ui);
}

inline bool operator!=(FLINT_Integer const& x, FLINT_Integer const& y) {
  return !(x == y);
}

inline bool operator!=(FLINT_Integer const& x, signed int const& y) {
  return !(x == y);
}

inline bool operator!=(FLINT_Integer const& x, signed long const& y) {
  return !(x == y);
}

inline bool operator!=(FLINT_Integer const& x, unsigned int const& y) {
  return !(x == y);
}

inline bool operator!=(FLINT_Integer const& x, unsigned long const& y) {
  return !(x == y);
}

inline bool operator<(FLINT_Integer const& x, FLINT_Integer const& y) {
  return (0 > fmpz_cmp(x.impl(), y.impl()));
}

inline bool operator<(FLINT_Integer const& x, signed int const& si) {
  return (0 > fmpz_cmp_si(x.impl(), si));
}

inline bool operator<(FLINT_Integer const& x, signed long const& si) {
  return (0 > fmpz_cmp_si(x.impl(), si));
}

inline bool operator<(FLINT_Integer const& x, unsigned int const& ui) {
  return (0 > fmpz_cmp_ui(x.impl(), ui));
}

inline bool operator<(FLINT_Integer const& x, unsigned long const& ui) {
  return (0 > fmpz_cmp_ui(x.impl(), ui));
}

inline bool operator>(FLINT_Integer const& x, FLINT_Integer const& y) {
  return (y < x);
}

inline bool operator>(FLINT_Integer const& x, signed int const& y) {
  return (y < x);
}

inline bool operator>(FLINT_Integer const& x, signed long const& y) {
  return (y < x);
}

inline bool operator>(FLINT_Integer const& x, unsigned int const& y) {
  return (y < x);
}

inline bool operator>(FLINT_Integer const& x, unsigned long const& y) {
  return (y < x);
}

inline bool operator>=(FLINT_Integer const& x, FLINT_Integer const& y) {
  return !(x < y);
}

inline bool operator>=(FLINT_Integer const& x, signed int const& y) {
  return !(x < y);
}

inline bool operator>=(FLINT_Integer const& x, signed long const& y) {
  return !(x < y);
}

inline bool operator>=(FLINT_Integer const& x, unsigned int const& y) {
  return !(x < y);
}

inline bool operator>=(FLINT_Integer const& x, unsigned long const& y) {
  return !(x < y);
}

inline bool operator<=(FLINT_Integer const& x, FLINT_Integer const& y) {
  return !(x > y);
}

inline bool operator<=(FLINT_Integer const& x, signed int const& y) {
  return !(x > y);
}

inline bool operator<=(FLINT_Integer const& x, signed long const& y) {
  return !(x > y);
}

inline bool operator<=(FLINT_Integer const& x, unsigned int const& y) {
  return !(x > y);
}

inline bool operator<=(FLINT_Integer const& x, unsigned long const& y) {
  return !(x > y);
}

inline FLINT_Integer& operator++(FLINT_Integer& x) {
  fmpz_add_ui(x.impl(), x.impl(), 1);
  return x;
}

inline FLINT_Integer& operator--(FLINT_Integer& x) {
  fmpz_sub_ui(x.impl(), x.impl(), 1);
  return x;
}

inline void abs_assign(FLINT_Integer& x) {
  fmpz_abs(x.impl(), x.impl());
}

inline FLINT_Integer abs(FLINT_Integer const& x) {
  FLINT_Integer res = x;
  abs_assign(res);
  return res;
}

inline void neg_assign(FLINT_Integer& x) {
  fmpz_neg(x.impl(), x.impl());
}

inline FLINT_Integer neg(FLINT_Integer const& x) {
  FLINT_Integer res = x;
  neg_assign(res);
  return res;
}

inline FLINT_Integer operator++(FLINT_Integer& x, int) {
  FLINT_Integer res(x);
  ++x;
  return res;
}

inline FLINT_Integer operator--(FLINT_Integer& x, int) {
  FLINT_Integer res(x);
  --x;
  return res;
}

inline void swap(FLINT_Integer& x, FLINT_Integer& y) noexcept {
  fmpz_swap(x.impl(), y.impl());
}

inline FLINT_Integer&
operator+=(FLINT_Integer& x, FLINT_Integer const& y) {
  fmpz_add(x.impl(), x.impl(), y.impl());
  return x;
}

inline FLINT_Integer&
operator+=(FLINT_Integer& x, signed int const& si) {
  x += FLINT_Integer(si);
  return x;
}

inline FLINT_Integer&
operator+=(FLINT_Integer& x, signed long const& si) {
  x += FLINT_Integer(si);
  return x;
}

inline FLINT_Integer&
operator+=(FLINT_Integer& x, unsigned int const& ui) {
  fmpz_add_ui(x.impl(), x.impl(), ui);
  return x;
}

inline FLINT_Integer&
operator+=(FLINT_Integer& x, unsigned long const& ui) {
  fmpz_add_ui(x.impl(), x.impl(), ui);
  return x;
}

inline FLINT_Integer&
operator-=(FLINT_Integer& x, FLINT_Integer const& y) {
  fmpz_sub(x.impl(), x.impl(), y.impl());
  return x;
}

inline FLINT_Integer&
operator-=(FLINT_Integer& x, signed int const& si) {
  x -= FLINT_Integer(si);
  return x;
}

inline FLINT_Integer&
operator-=(FLINT_Integer& x, signed long const& si) {
  x -= FLINT_Integer(si);
  return x;
}

inline FLINT_Integer&
operator-=(FLINT_Integer& x, unsigned int const& ui) {
  fmpz_sub_ui(x.impl(), x.impl(), ui);
  return x;
}

inline FLINT_Integer&
operator-=(FLINT_Integer& x, unsigned long const& ui) {
  fmpz_sub_ui(x.impl(), x.impl(), ui);
  return x;
}

inline FLINT_Integer&
operator*=(FLINT_Integer& x, FLINT_Integer const& y) {
  fmpz_mul(x.impl(), x.impl(), y.impl());
  return x;
}

inline FLINT_Integer&
operator*=(FLINT_Integer& x, signed int si) {
  fmpz_mul_si(x.impl(), x.impl(), si);
  return x;
}

inline FLINT_Integer&
operator*=(FLINT_Integer& x, signed long si) {
  fmpz_mul_si(x.impl(), x.impl(), si);
  return x;
}

inline FLINT_Integer&
operator*=(FLINT_Integer& x, unsigned int ui) {
  fmpz_mul_ui(x.impl(), x.impl(), ui);
  return x;
}

inline FLINT_Integer&
operator*=(FLINT_Integer& x, unsigned long ui) {
  fmpz_mul_ui(x.impl(), x.impl(), ui);
  return x;
}

inline FLINT_Integer&
operator/=(FLINT_Integer& x, FLINT_Integer const& y) {
  assert(y != 0);
  fmpz_tdiv_q(x.impl(), x.impl(), y.impl());
  return x;
}

inline FLINT_Integer&
operator/=(FLINT_Integer& x, signed int si) {
  assert(si != 0);
  fmpz_tdiv_q_si(x.impl(), x.impl(), si);
  return x;
}

inline FLINT_Integer&
operator/=(FLINT_Integer& x, signed long si) {
  assert(si != 0);
  fmpz_tdiv_q_si(x.impl(), x.impl(), si);
  return x;
}

inline FLINT_Integer&
operator/=(FLINT_Integer& x, unsigned int ui) {
  assert(ui != 0);
  fmpz_tdiv_q_ui(x.impl(), x.impl(), ui);
  return x;
}

inline FLINT_Integer&
operator/=(FLINT_Integer& x, unsigned long ui) {
  assert(ui != 0);
  fmpz_tdiv_q_ui(x.impl(), x.impl(), ui);
  return x;
}

inline FLINT_Integer&
operator%=(FLINT_Integer& x, FLINT_Integer const& y) {
  assert(y != 0);
  FLINT_Integer dummy;
  fmpz_tdiv_qr(dummy.impl(), x.impl(), x.impl(), y.impl());
  return x;
}

inline FLINT_Integer&
operator<<=(FLINT_Integer& x, unsigned long ui) {
  fmpz_mul_2exp(x.impl(), x.impl(), ui);
  return x;
}

inline FLINT_Integer&
operator>>=(FLINT_Integer& x, unsigned long ui) {
  fmpz_tdiv_q_2exp(x.impl(), x.impl(), ui);
  return x;
}

inline FLINT_Integer
operator+(FLINT_Integer const& x, FLINT_Integer const& y) {
  FLINT_Integer res;
  fmpz_add(res.impl(), x.impl(), y.impl());
  return res;
}

inline FLINT_Integer
operator+(FLINT_Integer const& x, signed int const& si) {
  return operator+(x, FLINT_Integer(si));
}

inline FLINT_Integer
operator+(FLINT_Integer const& x, signed long const& si) {
  return operator+(x, FLINT_Integer(si));
}

inline FLINT_Integer
operator+(FLINT_Integer const& x, unsigned int const& ui) {
  FLINT_Integer res;
  fmpz_add_ui(res.impl(), x.impl(), ui);
  return res;
}

inline FLINT_Integer
operator+(FLINT_Integer const& x, unsigned long const& ui) {
  FLINT_Integer res;
  fmpz_add_ui(res.impl(), x.impl(), ui);
  return res;
}

inline FLINT_Integer
operator-(FLINT_Integer const& x) {
  FLINT_Integer res = x;
  neg_assign(res);
  return res;
}

inline FLINT_Integer
operator-(FLINT_Integer const& x, FLINT_Integer const& y) {
  FLINT_Integer res;
  fmpz_sub(res.impl(), x.impl(), y.impl());
  return res;
}

inline FLINT_Integer
operator-(FLINT_Integer const& x, signed int const& si) {
  return operator-(x, FLINT_Integer(si));
}

inline FLINT_Integer
operator-(FLINT_Integer const& x, signed long const& si) {
  return operator-(x, FLINT_Integer(si));
}

inline FLINT_Integer
operator-(FLINT_Integer const& x, unsigned int const& ui) {
  FLINT_Integer res;
  fmpz_sub_ui(res.impl(), x.impl(), ui);
  return res;
}

inline FLINT_Integer
operator-(FLINT_Integer const& x, unsigned long const& ui) {
  FLINT_Integer res;
  fmpz_sub_ui(res.impl(), x.impl(), ui);
  return res;
}

inline FLINT_Integer
operator*(FLINT_Integer const& x, FLINT_Integer const& y) {
  FLINT_Integer res;
  fmpz_mul(res.impl(), x.impl(), y.impl());
  return res;
}

inline FLINT_Integer
operator*(FLINT_Integer const& x, signed int si) {
  FLINT_Integer res;
  fmpz_mul_si(res.impl(), x.impl(), si);
  return res;
}

inline FLINT_Integer
operator*(FLINT_Integer const& x, signed long si) {
  FLINT_Integer res;
  fmpz_mul_si(res.impl(), x.impl(), si);
  return res;
}

inline FLINT_Integer
operator*(FLINT_Integer const& x, unsigned int ui) {
  FLINT_Integer res;
  fmpz_mul_ui(res.impl(), x.impl(), ui);
  return res;
}

inline FLINT_Integer
operator*(FLINT_Integer const& x, unsigned long ui) {
  FLINT_Integer res;
  fmpz_mul_ui(res.impl(), x.impl(), ui);
  return res;
}

inline FLINT_Integer
operator/(FLINT_Integer const& x, FLINT_Integer const& y) {
  assert(y != 0);
  FLINT_Integer res;
  fmpz_tdiv_q(res.impl(), x.impl(), y.impl());
  return res;
}

inline FLINT_Integer
operator/(FLINT_Integer& x, signed int si) {
  assert(si != 0);
  FLINT_Integer res;
  fmpz_tdiv_q_si(res.impl(), x.impl(), si);
  return res;
}

inline FLINT_Integer
operator/(FLINT_Integer& x, signed long si) {
  assert(si != 0);
  FLINT_Integer res;
  fmpz_tdiv_q_si(res.impl(), x.impl(), si);
  return res;
}

inline FLINT_Integer
operator/(FLINT_Integer& x, unsigned int ui) {
  assert(ui != 0);
  FLINT_Integer res;
  fmpz_tdiv_q_ui(res.impl(), x.impl(), ui);
  return res;
}

inline FLINT_Integer
operator/(FLINT_Integer& x, unsigned long ui) {
  assert(ui != 0);
  FLINT_Integer res;
  fmpz_tdiv_q_ui(res.impl(), x.impl(), ui);
  return res;
}

inline FLINT_Integer
operator%(FLINT_Integer const& x, FLINT_Integer const& y) {
  FLINT_Integer res = x;
  res %= y;
  return res;
}

inline FLINT_Integer
operator<<(FLINT_Integer const& x, unsigned long ui) {
  FLINT_Integer res;
  fmpz_mul_2exp(res.impl(), x.impl(), ui);
  return res;
}

inline FLINT_Integer
operator>>(FLINT_Integer& x, unsigned int ui) {
  FLINT_Integer res;
  fmpz_tdiv_q_2exp(res.impl(), x.impl(), ui);
  return res;
}

inline int
compare(FLINT_Integer const& x, FLINT_Integer const& y) {
  return fmpz_cmp(x.impl(), y.impl());
}

inline void
add_mul_assign(FLINT_Integer& x,
               FLINT_Integer const& y, FLINT_Integer const& z) {
  fmpz_addmul(x.impl(), y.impl(), z.impl());
}

inline void
sub_mul_assign(FLINT_Integer& x,
               FLINT_Integer const& y, FLINT_Integer const& z) {
  fmpz_submul(x.impl(), y.impl(), z.impl());
}

inline void
gcd_assign(FLINT_Integer& x,
           FLINT_Integer const& y, FLINT_Integer const& z) {
  fmpz_gcd(x.impl(), y.impl(), z.impl());
}

inline void
lcm_assign(FLINT_Integer& x,
           FLINT_Integer const& y, FLINT_Integer const& z) {
  fmpz_lcm(x.impl(), y.impl(), z.impl());
}

inline void
exact_div_assign(FLINT_Integer& x,
                 FLINT_Integer const& y, FLINT_Integer const& z) {
  // Sets x to the quotient of y and z, assuming the division is exact.
  assert(z != 0);
  assert((y % z) == 0);
  fmpz_divexact(x.impl(), y.impl(), z.impl());
}

inline int
sgn(FLINT_Integer const& x) {
  return fmpz_sgn(x.impl());
}

inline void
mpz_set(mpz_t dst, const FLINT_Integer& src) {
  fmpz_get_mpz(dst, src.impl());
}

NOTHROW_DEFAULT_AND_MOVES(FLINT_Integer);

} // namespace pplite

#endif // !defined(pplite_FLINT_Integer_hh)
