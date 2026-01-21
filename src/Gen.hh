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

#ifndef pplite_Gen_hh
#define pplite_Gen_hh 1

#include "globals.hh"
#include "utils.hh"
#include "Bits.hh"
#include "Integer.hh"
#include "Linear_Expr.hh"
#include "Var.hh"

#include <algorithm>
#include <cassert>
#include <functional>
#include <iostream>
#include <tuple>
#include <utility>
#include <vector>

namespace pplite {

class Gen;
int compare(const Gen& x, const Gen& y);

class Gen {
public:
  enum Type {
    LINE,
    RAY,
    POINT,
    CLOSURE_POINT
  };

  struct Impl {
    Linear_Expr expr;
    Integer inhomo;
    Type type;
    Impl() noexcept : expr(), inhomo(1), type(POINT) {}
    Impl(Linear_Expr e, Integer d, Type t)
      : expr(std::move(e)), inhomo(std::move(d)), type(t) {}
  };

  Gen() = default;
  Gen(const Gen& g) = default;
  Gen(Gen&& g) = default;
  Gen& operator=(const Gen& g) = default;
  Gen& operator=(Gen&& g) = default;
  ~Gen() = default;

  Gen(Type t, Linear_Expr e, Integer d)
    : impl_(std::move(e), std::move(d), t) {
    if (impl().inhomo < 0) {
      neg_assign(impl().expr);
      neg_assign(impl().inhomo);
    }
    strong_normalize();
    assert(check_inv());
  }

  Impl& impl() { return impl_; }
  const Impl& impl() const { return impl_; }

  Type type() const { return impl().type; }
  void set_type(Type t) { impl().type = t; }

  bool is_line() const { return type() == LINE; }
  bool is_ray() const { return type() == RAY; }
  bool is_point() const { return type() == POINT; }
  bool is_closure_point() const { return type() == CLOSURE_POINT; }
  bool is_line_or_ray() const { return type() == LINE || type() == RAY; }
  bool is_point_or_closure_point() const { return !is_line_or_ray(); }

  Linear_Expr& linear_expr() { return impl().expr; }
  const Linear_Expr& linear_expr() const { return impl().expr; }

  const Integer& coeff(Var v) const { return impl().expr.get(v); }
  const Integer& divisor() const {
    assert(!is_line_or_ray());
    return impl().inhomo;
  }

  dim_type space_dim() const { return linear_expr().space_dim(); }
  void set_space_dim(dim_type dim) { linear_expr().set_space_dim(dim); }

  bool is_equal_to(const Gen& y) const {
    return type() == y.type()
      && linear_expr().is_equal_to(y.linear_expr())
      && impl().inhomo == y.impl().inhomo;
  }

  bool check_inv() const;

  void print(std::ostream& s) const;
  void print() const { print(std::cout); }

  void dump_type(std::ostream& s) const;
  void ascii_dump(std::ostream& s) const;

  bool load_type(std::istream& s);
  bool ascii_load(std::istream& s);

  void m_swap(Gen& y) noexcept {
    Impl& xi = impl();
    Impl& yi = y.impl();
    using std::swap;
    swap(xi.type, yi.type);
    swap(xi.expr, yi.expr);
    swap(xi.inhomo, yi.inhomo);
  }

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

  static constexpr Type LINE_OR_EQUALITY = LINE;
  bool is_line_or_equality() const { return is_line(); }

  void linear_combine(const Gen& y, dim_type dim) {
    impl().expr.linear_combine(y.linear_expr(), dim,
                               impl().inhomo, y.impl().inhomo);
    strong_normalize();
  }

  void sign_normalize() {
    if (is_line())
      impl().expr.sign_normalize();
  }

  void strong_normalize() {
    impl().expr.normalize(impl().inhomo);
    sign_normalize();
  }

  bool check_strong_normalized() const {
    Gen tmp = *this;
    tmp.strong_normalize();
    return compare(*this, tmp) == 0;
  }

private:
  Impl impl_;
}; // class Gen

NOTHROW_DEFAULT_AND_MOVES(Gen);

using Gens = std::vector<Gen>;

// lines, rays, closure_points, skel_points, nonskel_points
using Gens_Info = std::tuple<dim_type,dim_type,dim_type,dim_type,dim_type>;

inline Gen
line(Linear_Expr e) { return Gen(Gen::LINE, std::move(e), 0); }
inline Gen
ray(Linear_Expr e) { return Gen(Gen::RAY, std::move(e), 0); }
inline Gen
point(Linear_Expr e = Linear_Expr()) {
  return Gen(Gen::POINT, std::move(e), 1);
}
inline Gen
point(Linear_Expr e, Integer d) {
  return Gen(Gen::POINT, std::move(e), std::move(d));
}
inline Gen
closure_point(Linear_Expr e = Linear_Expr()) {
  return Gen(Gen::CLOSURE_POINT, std::move(e), 1);
}
inline Gen
closure_point(Linear_Expr e, Integer d) {
  return Gen(Gen::CLOSURE_POINT, std::move(e), std::move(d));
}

inline bool
operator==(const Gen& x, const Gen& y) { return x.is_equal_to(y); }
inline bool
operator!=(const Gen& x, const Gen& y) { return !(x == y); }

inline void swap(Gen& x, Gen& y) noexcept { x.m_swap(y); }

inline int
compare(const Gen& x, const Gen& y) {
  if (x.type() != y.type())
    return (x.type() < y.type()) ? -1 : 1;
  int res = compare(x.linear_expr(), y.linear_expr());
  return (res != 0)
    ? res
    : compare(x.impl().inhomo, y.impl().inhomo);
}

namespace IO_Operators {

std::ostream& operator<<(std::ostream& s, const Gen::Type& t);

inline std::ostream&
operator<<(std::ostream& s, const Gen& g) {
  g.print(s);
  return s;
}

std::ostream& operator<<(std::ostream& s, const Gens& gs);

} // namespace IO_Operators

inline bool
has_point(const Gens& gens) {
  return any_of(gens, std::mem_fn(&Gen::is_point));
}
inline bool
has_closure_point(const Gens& gens) {
  return any_of(gens, std::mem_fn(&Gen::is_closure_point));
}
inline bool
has_line_or_ray(const Gens& gens) {
  return any_of(gens, std::mem_fn(&Gen::is_line_or_ray));
}

template <typename Iter>
inline void
erase_space_dims(Gens& gs, Iter first, Iter last) {
  assert(first != last);
  for (auto& g : gs) {
    const auto dim = g.space_dim();
    Iter g_last = last;
    for ( ; first != g_last; --g_last) {
      Iter prev = g_last;
      --prev;
      if (*prev < dim)
        break;
    }
    erase_using_sorted_indices(g.linear_expr().impl(), first, g_last);
    g.strong_normalize();
  }
}

inline Index_Set
invalid_lines(const Gens& gs) {
  Index_Set res;
  for (auto i : bwd_index_range(gs)) {
    if (gs[i].is_line() && gs[i].linear_expr().is_zero())
      res.set(i);
  }
  return res;
}

inline Index_Set
invalid_rays(const Gens& gs) {
  Index_Set res;
  for (auto i : bwd_index_range(gs)) {
    if (gs[i].is_ray() && gs[i].linear_expr().is_zero())
      res.set(i);
  }
  return res;
}

namespace detail {

inline void
erase_higher_dims(Gens& gs, dim_type d) {
  for (auto& g : gs) {
    if (g.space_dim() > d) {
      g.set_space_dim(d);
      g.strong_normalize();
    }
  }
}

void add_as_rays(Gens gens, Gens& rays);

template <typename Indices>
Gen
materialize(const Indices& is, const Gens& gs) {
  Linear_Expr expr;
  Integer div;
  for (auto i : is) {
    const auto& g = gs[i];
    expr += g.linear_expr();
    if (g.is_closure_point())
      div += g.divisor();
  }
  return point(std::move(expr), std::move(div));
}

} // namespace detail

} // namespace pplite

#endif // !defined(pplite_Gen_hh)
