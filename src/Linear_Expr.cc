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

#include "pplite-config.h"
#include "Linear_Expr.hh"
#include "Low_Level_Stats.hh"
#include "ascii_dump_load.hh"

#include <iostream>
#include <algorithm>

namespace pplite {

int
compare(Linear_Expr const& x, Linear_Expr const& y) {
  const dim_type min_dim = std::min(x.space_dim(), y.space_dim());
  for (auto i : range(min_dim)) {
    auto res = compare(x.get(i), y.get(i));
    if (res != 0)
      return res;
  }
  for (auto i : range(min_dim, x.space_dim())) {
    auto res = sgn(x.get(i));
    if (res != 0)
      return res;
  }
  for (auto i : range(min_dim, y.space_dim())) {
    auto res = -sgn(y.get(i));
    if (res != 0)
      return res;
  }
  return 0;
}

Integer
Linear_Expr::gcd(dim_type first, dim_type last) const {
  assert(0 <= first && first <= last && last <= space_dim());
  auto fnz = first_nonzero(first, last);
  if (fnz == last)
    return Integer::zero();
  const auto& x = impl();
  Integer res = x[fnz];
  if (sgn(res) < 0)
    neg_assign(res);
  if (res == 1)
    return res;
  for (auto i : range(fnz+1, last)) {
    if (x[i] != 0) {
      gcd_assign(res, res, x[i]);
      if (res == 1)
        return res;
    }
  }
  return res;
}

void
Linear_Expr::normalize(Integer& inhomo) {
  Integer gcd_in = gcd(0, space_dim());
  if (inhomo != 0 && gcd_in != 1) {
    if (gcd_in != 0)
      gcd_assign(gcd_in, gcd_in, inhomo);
    else {
      gcd_in = inhomo;
      if (sgn(gcd_in) < 0)
        neg_assign(gcd_in);
    }
  }
  if (gcd_in == 0 || gcd_in == 1)
    return;

  // Divide the coefficients by the GCD.
  exact_div_assign(inhomo, inhomo, gcd_in);
  auto& x = impl();
  for (auto i : bwd_index_range(x)) {
    auto& x_i = x[i];
    exact_div_assign(x_i, x_i, gcd_in);
  }
}

void
Linear_Expr::combine(Linear_Expr& x, Integer& x_inhomo,
                     const Linear_Expr& y, const Integer& y_inhomo,
                     const Integer& x_coeff, const Integer& y_coeff) {
#if PPLITE_LOW_LEVEL_COUNTERS
  LLOp_Stats::bump_linear_comb();
#endif
  assert(x_coeff != 0 && y_coeff != 0);
  static PPLITE_TLS Integer x_coprime;
  static PPLITE_TLS Integer y_coprime;
  get_coprimes(x_coeff, y_coeff, x_coprime, y_coprime);
  neg_assign(x_coprime);

  const auto x_sd = x.space_dim();
  const auto y_sd = y.space_dim();
  const auto min_dim = std::min(x_sd, y_sd);
  if (x_sd < y_sd)
    x.set_space_dim(y_sd);
  // x_sd still refers to old space dim (meant!)

  auto combine_p1p1 = [&] () {
    x_inhomo += y_inhomo;
    for (auto i : range(min_dim))
      x[i] += y[i];
    if (y_sd > min_dim) {
      for (auto i : range(min_dim, y_sd))
        x[i] = y[i];
    }
  }; // combine_p1p1

  auto combine_p1m1 = [&] () {
    neg_assign(x_inhomo);
    x_inhomo += y_inhomo;
    for (auto i : range(min_dim)) {
      neg_assign(x[i]);
      x[i] += y[i];
    }
    if (x_sd > min_dim) {
      for (auto i : range(min_dim, x_sd))
        neg_assign(x[i]);
    } else if (y_sd > min_dim) {
      for (auto i : range(min_dim, y_sd))
        x[i] = y[i];
    }
  }; // combine_p1m1

  auto combine_m1p1 = [&] () {
    x_inhomo -= y_inhomo;
    for (auto i : range(min_dim))
      x[i] -= y[i];
    if (y_sd > min_dim) {
      for (auto i : range(min_dim, y_sd))
        x[i] -= y[i];
    }
  }; // combine_m1p1

  auto combine_m1m1 = [&] () {
    neg_assign(x_inhomo);
    x_inhomo -= y_inhomo;
    for (auto i : range(min_dim)) {
      neg_assign(x[i]);
      x[i] -= y[i];
    }
    if (x_sd > min_dim) {
      for (auto i : range(min_dim, x_sd))
        neg_assign(x[i]);
    } else if (y_sd > min_dim) {
      for (auto i : range(min_dim, y_sd))
        x[i] -= y[i];
    }
  }; // combine_m1m1

  auto combine_Zp1 = [&] () {
    add_mul_assign(x_inhomo, y_inhomo, x_coprime);
    for (auto i : range(min_dim))
      add_mul_assign(x[i], y[i], x_coprime);
    if (y_sd > min_dim) {
      for (auto i : range(min_dim, y_sd)) {
        x[i] = y[i];
        x[i] *= x_coprime;
      }
    }
  }; // combine_Zp1

  auto combine_Zm1 = [&] () {
    neg_assign(x_inhomo);
    add_mul_assign(x_inhomo, y_inhomo, x_coprime);
    for (auto i : range(min_dim)) {
      neg_assign(x[i]);
      add_mul_assign(x[i], y[i], x_coprime);
    }
    if (x_sd > min_dim) {
      for (auto i : range(min_dim, x_sd))
        neg_assign(x[i]);
    } else if (y_sd > min_dim) {
      for (auto i : range(min_dim, y_sd)) {
        x[i] = y[i];
        x[i] *= x_coprime;
      }
    }
  }; // combine_Zm1

  auto combine_p1Z = [&] () {
    x_inhomo *= y_coprime;
    x_inhomo += y_inhomo;
    for (auto i : range(min_dim)) {
      x[i] *= y_coprime;
      x[i] += y[i];
    }
    if (x_sd > min_dim) {
      for (auto i : range(min_dim, x_sd))
        x[i] *= y_coprime;
    } else if (y_sd > min_dim) {
      for (auto i : range(min_dim, y_sd))
        x[i] = y[i];
    }
  }; // combine_p1Z

  auto combine_m1Z = [&] () {
    x_inhomo *= y_coprime;
    x_inhomo -= y_inhomo;
    for (auto i : range(min_dim)) {
      x[i] *= y_coprime;
      x[i] -= y[i];
    }
    if (x_sd > min_dim) {
      for (auto i : range(min_dim, x_sd))
        x[i] *= y_coprime;
    } else if (y_sd > min_dim) {
      for (auto i : range(min_dim, y_sd))
        x[i] -= y[i];
    }
  }; // combine_p1Z

  auto combine_ZZ = [&] () {
    x_inhomo *= y_coprime;
    add_mul_assign(x_inhomo, y_inhomo, x_coprime);
    for (auto i : range(min_dim)) {
      x[i] *= y_coprime;
      add_mul_assign(x[i], y[i], x_coprime);
    }
    if (x_sd > min_dim) {
      for (auto i : range(min_dim, x_sd))
        x[i] *= y_coprime;
    } else if (y_sd > min_dim) {
      for (auto i : range(min_dim, y_sd)) {
        x[i] = y[i];
        x[i] *= x_coprime;
      }
    }
  }; // combine_ZZ

  // Silly way to use a switch rather than a chain of if-then-else.
  const int mask
    = (x_coprime.is_pm1() << 3)
    + (y_coprime.is_pm1() << 2)
    + (x_coprime.is_one() << 1)
    + (y_coprime.is_one() << 0);
  switch (mask) {
  case 15:
    combine_p1p1();
    break;
  case 14:
    combine_p1m1();
    break;
  case 13:
    combine_m1p1();
    break;
  case 12:
    combine_m1m1();
    break;
  case 10:
    combine_p1Z();
    break;
  case 8:
    combine_m1Z();
    break;
  case 5:
    combine_Zp1();
    break;
  case 4:
    combine_Zm1();
    break;
  case 0:
    combine_ZZ();
    break;
  default:
    PPLITE_UNREACH;
  }
}

void Linear_Expr::print(std::ostream& os) const {
  if (is_zero()) {
    os << "0";
    return;
  }

  bool first = true;
  for (auto i : dim_range(*this)) {
    auto s = sgn(row[i]);
    if (s == 0)
      continue;
    if (s > 0) {
      if (!first)
        os << " + ";
      if (row[i] != 1) {
        row[i].print(os);
        os << "*";
      }
    } else {
      assert(s < 0);
      if (first)
        os << "-";
      else
        os << " - ";
      if (row[i] != -1) {
        neg(row[i]).print(os);
        os << "*";
      }
    }
    first = false;
    Var(i).print(os);
  }
}

void
Linear_Expr::ascii_dump(std::ostream& s) const {
  using IO_Operators::operator<<;
  s << "dim " << space_dim() << " : ";
  for (auto i : dim_range(*this)) {
    get(i).print(s);
    s << ' ';
  }
}

bool
Linear_Expr::ascii_load(std::istream& s) {
  dim_type dim = 0;
  if (!(ascii_load_string(s, "dim")
        && (s >> dim)
        && (dim >= 0)
        && ascii_load_string(s, ":")))
    return false;
  row.resize(dim);
  for (auto& i : row) {
    if (!i.ascii_load(s))
      return false;
  }
  return true;
}

} // namespace pplite
