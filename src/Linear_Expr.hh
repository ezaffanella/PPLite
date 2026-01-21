/* PPLite: a lightweight library for convex polyhedra derived from PPL.
   Copyright (C) 2001-2010 Roberto Bagnara <bagnara@cs.unipr.it>
   Copyright (C) 2010-2018 BUGSENG srl (http://bugseng.com)
   Copyright (C) 2018-2026 Enea Zaffanella <enea.zaffanella@unipr.it>

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

#ifndef pplite_Linear_Expr_hh
#define pplite_Linear_Expr_hh 1

#include "globals.hh"
#include "Bits.hh"
#include "Integer.hh"
#include "Var.hh"
#include "utils.hh"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <utility>
#include <vector>

namespace pplite {

class Linear_Expr {
private:
  using Impl = std::vector<Integer>;
  Impl row;

public:
  explicit Linear_Expr(dim_type dim) : row(dim) {}

  // Not explicit: should behave as implicit conversion.
  Linear_Expr(Var v) : row(v.space_dim()) { row[v.id()] = 1; }
  Linear_Expr(Linear_Expr const& e, dim_type dim)
    : row(dim) {
    assert(e.space_dim() <= dim);
    std::copy(e.row.cbegin(), e.row.cend(), row.begin());
  }

  Impl& impl() { return row; }
  Impl const& impl() const { return row; }

  void print(std::ostream& os) const;
  void print() const { print(std::cout); }

  Linear_Expr() = default;
  Linear_Expr(Linear_Expr const& e) = default;
  Linear_Expr(Linear_Expr&& e) = default;
  Linear_Expr& operator=(Linear_Expr const& e) = default;
  Linear_Expr& operator=(Linear_Expr&& e) = default;
  ~Linear_Expr() = default;

  using iterator = Impl::iterator;
  iterator begin() { return impl().begin(); }
  iterator end() { return impl().end(); }

  using const_iterator = Impl::const_iterator;
  const_iterator begin() const { return impl().begin(); }
  const_iterator end() const { return impl().end(); }
  const_iterator cbegin() const { return impl().cbegin(); }
  const_iterator cend() const { return impl().cend(); }

  dim_type space_dim() const { return row.size(); }
  void set_space_dim(dim_type dim) { row.resize(dim); }

  void swap_space_dims(dim_type i, dim_type j) {
    const auto min_dim = 1 + std::max(i, j);
    if (space_dim() < min_dim)
      set_space_dim(min_dim);
    using std::swap;
    swap(row[i], row[j]);
  }

  // Note: [first, last) has to be an ordered set of space dims.
  template <typename Iter>
  void remove_space_dims(Iter first, Iter last) {
    erase_using_sorted_indices(row, first, last);
  }

  // Shift by n positions the coefficients of variables,
  // starting from space dim `start'.
  void shift_space_dims(dim_type start, dim_type n) {
    dim_type position = space_dim() - 1;
    set_space_dim(space_dim() + n);
    while (position >= start) {
      using std::swap;
      swap(row[position], row[position + n]);
      assert(row[position] == 0);
      --position;
    }
  }

  void shift_space_dims(Var start, dim_type n) {
    shift_space_dims(start.id(), n);
  }

  void permute_space_dims_cycle(const Dims& cycle, dim_type d) {
    if (cycle.size() <= 1)
      return;
    if (cycle.size() == 2) {
      swap_space_dims(cycle[0], cycle[1]);
      return;
    }
    if (space_dim() < d)
      set_space_dim(d);
    Integer tmp = row[cycle.back()];
    using std::swap;
    for (auto i : cycle)
      swap(row[i], tmp);
  }

  bool is_zero() const {
    return all_zeroes(0, space_dim());
  }

  Index_Set non_zeroes() const {
    Index_Set res;
    for (auto i : bwd_dim_range(*this)) {
      if (not row[i].is_zero())
        res.set(i);
    }
    return res;
  }

  std::pair<Index_Set, Index_Set> pos_neg() const {
    Index_Set pos, neg;
    for (auto i : bwd_dim_range(*this)) {
      auto s = sgn(row[i]);
      if (s == 0)
        continue;
      else if (s == 1)
        pos.set(i);
      else
        neg.set(i);
    }
    return { pos, neg };
  }

  void ascii_dump(std::ostream& s) const;
  bool ascii_load(std::istream& s);

  void m_swap(Linear_Expr& y) noexcept {
    using std::swap;
    swap(row, y.row);
  }

  void normalize(Integer& inhomo);

  void sign_normalize() {
    auto fnz = first_nonzero();
    if (fnz < space_dim() && row[fnz] < 0)
      negate(fnz, space_dim());
  }

  void sign_normalize(Integer& inhomo) {
    auto fnz = first_nonzero();
    if (fnz < space_dim() && row[fnz] < 0) {
      negate(fnz, space_dim());
      neg_assign(inhomo);
    }
  }

  Integer& operator[](dim_type dim) {
    assert(0 <= dim && dim < space_dim());
    return row[dim];
  }
  const Integer& operator[](dim_type dim) const {
    assert(0 <= dim && dim < space_dim());
    return row[dim];
  }

  Integer const& get(dim_type dim) const {
    return dim < space_dim() ? row[dim] : Integer::zero();
  }
  Integer const& get(Var v) const { return get(v.id()); }

  void set(dim_type dim, Integer const& n) {
    if (dim + 1 > space_dim()) row.resize(dim + 1);
    row[dim] = n;
  }
  void set(Var v, Integer const& n) { set(v.id(), n); }

  // bool all_zeroes(Vars_Set const& vars) const;
  bool all_zeroes(dim_type first, dim_type last) const {
    assert(0 <= first && first <= last && last <= space_dim());
    return std::all_of(row.cbegin() + first, row.cbegin() + last,
                       [](Integer const& i) { return i.is_zero(); });
  }

  dim_type num_zeroes(dim_type first, dim_type last) const {
    assert(0 <= first && first <= last && last <= space_dim());
    return std::count_if(row.cbegin() + first, row.cbegin() + last,
                         [](const Integer& i) { return i.is_zero(); });
  }

  Integer gcd(dim_type first, dim_type last) const;

  static void
  combine(Linear_Expr& x, Integer& x_inhomo,
          const Linear_Expr& y, const Integer& y_inhomo,
          const Integer& x_coeff, const Integer& y_coeff);

  void linear_combine(const Linear_Expr& y, dim_type dim,
                      Integer& x_inhomo, const Integer& y_inhomo) {
    auto& x = *this;
    assert(dim < x.space_dim() && dim < y.space_dim());
    combine(x, x_inhomo, y, y_inhomo, x[dim], y[dim]);
  }

  void mul_assign(Integer const& n, dim_type first, dim_type last) {
    assert(first >= 0 && first <= last && last <= space_dim());
    std::for_each(row.begin() + first, row.begin() + last,
                  [&n](Integer& i) { i *= n; });
  }

  dim_type last_nonzero() const {
    return last_nonzero(0, space_dim());
  }

  // Returns the index of the last nonzero element in [first,last),
  // or last if there are no nonzero elements.
  dim_type last_nonzero(dim_type first, dim_type last) const {
    assert(0 <= first && first <= last && last <= space_dim());
    auto i = last;
    while (i != first) {
      --i;
      if (row[i] != 0)
        return i;
    }
    return last;
  }

  // Returns the index of the first nonzero element in [first,last),
  // or last if there are no nonzero elements.
  dim_type first_nonzero(dim_type first, dim_type last) const {
    assert(0 <= first && first <= last && last <= space_dim());
    for (auto i : range(first, last)) {
      if (row[i] != 0)
        return i;
    }
    return last;
  }

  dim_type first_nonzero() const {
    return first_nonzero(0, space_dim());
  }

  bool all_zeroes_except(Var var) const {
    return get(var.id()) != 0
      && all_zeroes(0, var.id())
      && all_zeroes(var.id() + 1, space_dim());
  }

  static bool
  is_equal_to(const Linear_Expr& x, const Linear_Expr& y,
              dim_type first, dim_type last) {
    assert(0 <= first && first <= last);
    assert(last <= std::min(x.space_dim(), y.space_dim()));
    return std::equal(x.cbegin() + first, x.cbegin() + last,
                      y.cbegin() + first);
  }

  bool is_equal_to(const Linear_Expr& y) const {
    const auto& x = *this;
    auto min_dim = std::min(x.space_dim(), y.space_dim());
    return is_equal_to(x, y, 0, min_dim)
      && x.all_zeroes(min_dim, x.space_dim())
      && y.all_zeroes(min_dim, y.space_dim());
  }

  static bool
  is_equal_to(const Linear_Expr& x, const Linear_Expr& y,
              const Integer& x_coeff, const Integer& y_coeff,
              dim_type first, dim_type last) {
    assert(0 <= first && first <= last);
    assert(last <= std::min(x.space_dim(), y.space_dim()));
    assert(x_coeff != 0 && y_coeff != 0);
    return std::equal(x.cbegin() + first, x.cbegin() + last,
                      y.cbegin() + first,
                      [&x_coeff, &y_coeff](const Integer& xi,
                                           const Integer& yi) {
                        return xi * x_coeff == yi * y_coeff;
                      });
  }

  bool is_equal_to(const Linear_Expr& y,
                   const Integer& coeff, const Integer& y_coeff) const {
    if (y_coeff == 0)
      return (coeff == 0) || is_zero();
    if (coeff == 0)
      return y.is_zero();
    const auto& x = *this;
    auto min_dim = std::min(x.space_dim(), y.space_dim());
    return is_equal_to(x, y, coeff, y_coeff, 0, min_dim)
      && x.all_zeroes(min_dim, x.space_dim())
      && y.all_zeroes(min_dim, y.space_dim());
  }

  void negate(dim_type first, dim_type last) {
    assert(0 <= first && first <= last && last <= space_dim());
    std::for_each(row.begin() + first, row.begin() + last,
                  [](Integer& i) { neg_assign(i); });
  }

}; // class Linear_Expr

inline Linear_Expr&
operator+=(Linear_Expr& e1, Linear_Expr const& e2) {
  if (e1.space_dim() < e2.space_dim())
    e1.set_space_dim(e2.space_dim());
  for (auto i : bwd_dim_range(e2))
    e1[i] += e2[i];
  return e1;
}

inline Linear_Expr&
operator+=(Linear_Expr& e, Var v) {
  if (e.space_dim() < v.space_dim())
    e.set_space_dim(v.space_dim());
  ++e[v.id()];
  return e;
}

inline Linear_Expr&
operator-=(Linear_Expr& e1, Linear_Expr const& e2) {
  if (e1.space_dim() < e2.space_dim())
    e1.set_space_dim(e2.space_dim());
  for (auto i : bwd_dim_range(e2))
    e1[i] -= e2[i];
  return e1;
}

inline Linear_Expr&
operator-=(Linear_Expr& e, Var v) {
  if (e.space_dim() < v.space_dim())
    e.set_space_dim(v.space_dim());
  --e[v.id()];
  return e;
}

inline Linear_Expr&
operator*=(Linear_Expr& e, Integer const& n) {
  for (auto& i : e)
    i *= n;
  return e;
}

inline void
neg_assign(Linear_Expr& e) {
  for (auto& i : e)
    neg_assign(i);
}

inline void
add_mul_assign(Linear_Expr& e, Integer const& n, Var v) {
  if (v.space_dim() > e.space_dim())
    e.set_space_dim(v.space_dim());
  e[v.id()] += n;
}

inline void
sub_mul_assign(Linear_Expr& e, Integer const& n, Var v) {
  if (e.space_dim() < v.space_dim())
    e.set_space_dim(v.space_dim());
  e[v.id()] -= n;
}

inline void
add_mul_assign(Linear_Expr& e1, Integer const& n, Linear_Expr const& e2) {
  if (e2.space_dim() > e1.space_dim())
    e1.set_space_dim(e2.space_dim());
  for (auto i : bwd_dim_range(e2))
    add_mul_assign(e1[i], n, e2[i]);
}

inline void
sub_mul_assign(Linear_Expr& e1, Integer const& n, Linear_Expr const& e2) {
  if (e2.space_dim() > e1.space_dim())
    e1.set_space_dim(e2.space_dim());
  for (auto i : bwd_dim_range(e2))
    sub_mul_assign(e1[i], n, e2[i]);
}

inline Linear_Expr
operator+(Linear_Expr e1, Linear_Expr const& e2) { return e1 += e2; }

inline Linear_Expr
operator+(Var v, Var w) {
  Var v_max(std::max(v.id(), w.id()));
  Var v_min(std::min(v.id(), w.id()));
  Linear_Expr e(v_max);
  return e += v_min;
}

inline Linear_Expr
operator+(Linear_Expr e, Var v) { return e += v; }

inline Linear_Expr
operator+(Var v, Linear_Expr e) { return e += v; }

inline Linear_Expr
operator+(Linear_Expr e) { return e; }

inline Linear_Expr
operator-(Linear_Expr e) {
  neg_assign(e);
  return e;
}

inline Linear_Expr
operator-(Linear_Expr e1, Linear_Expr const& e2) { return e1 -= e2; }

inline Linear_Expr
operator-(Var v, Var w) {
  auto max_sd = std::max(v.id(), w.id());
  Linear_Expr e = Var(max_sd);
  if (max_sd == v.id())
    e -= w;
  else {
    neg_assign(e);
    e += v;
  }
  return e;
}

inline Linear_Expr
operator-(Var v, Linear_Expr e) {
  neg_assign(e);
  return e += v;
}

inline Linear_Expr
operator-(Linear_Expr e, Var v) { return e -= v; }

inline Linear_Expr
operator*(Integer const& n, Linear_Expr e) { return e *= n; }

inline Linear_Expr
operator*(Linear_Expr e, Integer const& n) { return n * std::move(e); }

int compare(Linear_Expr const& x, Linear_Expr const& y);

inline bool
is_interval_expr(const Linear_Expr& e) {
  auto i = e.first_nonzero();
  return (i == e.space_dim()) || e.all_zeroes(i+1, e.space_dim());
}

inline bool
is_proper_interval_expr(const Linear_Expr& e) {
  auto i = e.first_nonzero();
  return (i < e.space_dim()) && e.all_zeroes(i+1, e.space_dim());
}

namespace IO_Operators {

inline std::ostream&
operator<<(std::ostream& s, Linear_Expr const& e) {
  e.print(s);
  return s;
}

} // namespace IO_Operators

inline void swap(Linear_Expr& x, Linear_Expr& y) noexcept {
  x.m_swap(y);
}

NOTHROW_DEFAULT_AND_MOVES(Linear_Expr);

using Linear_Exprs = std::vector<Linear_Expr>;

} // namespace pplite

#endif // !defined(pplite_Linear_Expr_hh)
