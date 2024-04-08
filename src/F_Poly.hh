/* PPLite: a lightweight library for convex polyhedra derived from PPL.
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

#ifndef pplite_F_Poly_hh
#define pplite_F_Poly_hh 1

#include "globals.hh"
#include "Poly.hh"
#include "Poly_min.hh"
#include "Poly_templ.hh"

#include <cassert>
#include <functional>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <numeric>
#include <vector>

namespace pplite {

struct F_Poly_Impl {
  using Block = std::vector<dim_type>;
  using Blocks = std::vector<Block>;
  using Factor = Poly;
  using Factors = std::vector<Factor>;

  dim_type dim;
  Topol topol;
  bool empty;
  bool is_normalized;
  Itvs itvs;
  Blocks blocks;
  Factors factors;

  F_Poly_Impl(dim_type d, Spec_Elem s, Topol t)
    : dim(d),
      topol(t),
      empty(s == Spec_Elem::EMPTY),
      is_normalized(true),
      itvs(),
      blocks(),
      factors() {
    if (empty)
      return;
    itvs.resize(dim);
  }

  // Note: in itvs, those dimensions occurring in blocks are mapped
  // to empty intervals; use this proxy to iterate on all
  // *proper* intervals, i.e., intervals of non-block dimensions.
  struct Proper_Itvs {
    const Itvs* ptr;

    explicit Proper_Itvs(const Itvs& itvs) : ptr(&itvs) {}

    using value_type = Itvs::value_type;
    const value_type* get_value_ptr(dim_type pos) const {
      return &((*ptr)[pos]);
    }
    bool is_skippable(dim_type pos) const {
      return get_value_ptr(pos)->is_empty();
    }
    dim_type end_pos() const { return num_rows(*ptr); }
    dim_type size() const { return std::distance(begin(), end()); }

    using const_iterator = Proxy_Iter<Proper_Itvs>;
    const_iterator begin() const { return const_iterator(this, false); }
    const_iterator end() const { return const_iterator(this, true); }
  }; // Proper_Itvs

  Proper_Itvs proper_itvs() const { return Proper_Itvs(itvs); }
}; // F_Poly_Impl

class F_Poly : private F_Poly_Impl {
public:
  using Impl = F_Poly_Impl;
  Impl& impl() { return *this; }
  const Impl& impl() const { return *this; }

  explicit
  F_Poly(dim_type d = 0, Spec_Elem s = Spec_Elem::UNIVERSE,
         Topol t = get_default_topology())
    : F_Poly_Impl(d, s, t) {}
  F_Poly(dim_type d, Topol t, Spec_Elem s = Spec_Elem::UNIVERSE)
    : F_Poly(d, s, t) {}
  F_Poly(Spec_Elem s, dim_type d, Topol t = get_default_topology())
    : F_Poly(d, s, t) {}
  F_Poly(Spec_Elem s, Topol t, dim_type d)
    : F_Poly(d, s, t) {}
  F_Poly(Topol t, dim_type d, Spec_Elem s = Spec_Elem::UNIVERSE)
    : F_Poly(d, s, t) {}
  F_Poly(Topol t, Spec_Elem s, dim_type d)
    : F_Poly(d, s, t) {}

  F_Poly(const F_Poly& y) = default;
  F_Poly& operator=(const F_Poly& y) = default;
  F_Poly(F_Poly&& y) = default;
  F_Poly& operator=(F_Poly&& y) = default;
  ~F_Poly() = default;

  /* Predicates */
  bool check_inv() const;
  bool is_necessarily_closed() const { return topology() == Topol::CLOSED; }
  bool is_empty() const { minimize(); return empty; }
  bool is_minimized() const {
    return empty || all_of(factors, std::mem_fn(&Factor::is_minimized));
  }
  bool is_universe() const {
    return !empty
      && all_of(proper_itvs(), std::mem_fn(&Itv::is_universe))
      && all_of(factors, std::mem_fn(&Factor::is_universe));
  }
  bool is_topologically_closed() const {
    return is_necessarily_closed()
      || empty
      // TODO: if Itv will support NNC, check them here.
      || all_of(factors, std::mem_fn(&Factor::is_topologically_closed));
  };
  bool is_bounded() const {
    return empty
      || (all_of(proper_itvs(), std::mem_fn(&Itv::is_bounded))
          &&
          all_of(factors, std::mem_fn(&Factor::is_bounded)));
  }
  bool is_bounded_expr(bool from_below, const Linear_Expr& expr) const;
  bool constrains(Var var) const;
  bool equals(const F_Poly& y) const;
  bool contains(const F_Poly& y) const;
  bool strictly_contains(const F_Poly& y) const;
  bool is_disjoint_from(const F_Poly& y) const;
  bool boxed_contains(const F_Poly& y) const;

  /* Queries */
  Topol topology() const { return topol; }
  dim_type space_dim() const { return dim; }
  dim_type affine_dim() const;
  Poly_Con_Rel relation_with(const Con& c) const;
  Poly_Gen_Rel relation_with(const Gen& g) const;
  bool min(const Affine_Expr& ae, Rational& value,
           bool* included_ptr = nullptr, Gen* g_ptr = nullptr) const;
  bool max(const Affine_Expr& ae, Rational& value,
           bool* included_ptr = nullptr, Gen* g_ptr = nullptr) const {
    bool res = min(-ae, value, included_ptr, g_ptr);
    neg_assign(value);
    return res;
  }
  Itv get_bounds(Var var) const;
  Itv get_bounds(const Affine_Expr& ae) const;
  Itv get_bounds(const Itv_Expr& ie) const;
  Index_Set get_unconstrained() const;

  size_t hash() const;

  void print(std::ostream& os) const;
  void print() const { print(std::cout); }

  using Cons_Proxy = Cont_Proxy<Cons>;
  Cons_Proxy cons() const { return Cons_Proxy(copy_cons()); }
  using Gens_Proxy = Cont_Proxy<Gens>;
  Gens_Proxy gens() const { return Gens_Proxy(copy_gens()); }
  Cons copy_cons() const;
  Gens copy_gens() const;
  Cons_Proxy normalized_cons() const { normalize(); return cons(); }
  Cons copy_normalized_cons() const  { normalize(); return copy_cons(); }
  Gens_Info gens_info() const;
  dim_type num_min_cons() const;
  dim_type num_min_gens() const;

  BBox get_bounding_box() const;

  void collapse(dim_type) { /* nothing to do */ }
  dim_type num_disjuncts() const { return is_empty() ? 0 : 1; }
  Cons_Proxy disjunct_cons(dim_type n) const {
    (void) n;
    assert(n == 0);
    return cons();
  }
  bool geom_covers(const F_Poly& y) const { return contains(y); }

  /* Modifiers */
  void m_swap(F_Poly& y) noexcept;
  void set_empty();
  void set_universe();
  void set_topology(Topol t);

  void add_con(const Con& c);
  void add_cons(const Cons& cs) { add_cons(cs.begin(), cs.end()); }
  // FIXME: reconsider if normalization should be the default.
  template <typename Iter>
  void add_cons(Iter first, Iter last, bool do_normalize = true);

  void add_gen(const Gen& g) { add_gens(&g, &g + 1); }
  void add_gens(const Gens& gs) { add_gens(gs.begin(), gs.end()); }
  template <typename Iter>
  void add_gens(Iter first, Iter last);

  void topological_closure_assign();

  void unconstrain(const Index_Set& vars);
  void unconstrain(Var var) {
    unconstrain(Index_Set(var.id()));
  }
  template <typename Iter>
  void unconstrain(Iter first, Iter last) {
    unconstrain(Index_Set(first, last));
  }

  void intersection_assign(const F_Poly& y);
  void join_assign(const F_Poly& y) { poly_hull_assign(y); }
  void poly_hull_assign(const F_Poly& y);
  void con_hull_assign(const F_Poly& y, bool boxed = false) {
    con_hull(*this, &y, &y + 1, boxed);
  }
  void poly_difference_assign(const F_Poly& y) {
    auto& x = *this;
    x = detail::poly_difference(x, y);
  }

  void affine_image(Var var, const Linear_Expr& expr,
                    const Integer& inhomo = Integer::zero(),
                    const Integer& den = Integer::one());
  void affine_preimage(Var var, const Linear_Expr& expr,
                       const Integer& inhomo = Integer::zero(),
                       const Integer& den = Integer::one());
  void
  parallel_affine_image(const Vars& vars,
                        const Linear_Exprs& exprs,
                        const Integers& inhomos, const Integers& dens);

  void widening_assign(const F_Poly& y, const Cons* upto_ptr,
                       Widen_Impl w_impl, Widen_Spec w_spec);
  void widening_assign(const F_Poly& y,
                       Widen_Impl w_impl = get_widen_impl(),
                       Widen_Spec w_spec = get_widen_spec()) {
    widening_assign(y, nullptr, w_impl, w_spec);
  }
  void widening_assign(const F_Poly& y, const Cons& upto_cons,
                       Widen_Impl w_impl = get_widen_impl(),
                       Widen_Spec w_spec = get_widen_spec()) {
    widening_assign(y, &upto_cons, w_impl, w_spec);
  }

  void time_elapse_assign(const F_Poly&);

private:
  F_Poly split_aux(const Con& c, Topol t, bool integral);
public:
  F_Poly split(const Con& c, Topol t) { return split_aux(c, t, false); }
  F_Poly split(const Con& c) { return split(c, topology()); }
  F_Poly integral_split(const Con& c) {
    return split_aux(c, topology(), true);
  }

  /* Change of space dim */
  void add_space_dims(dim_type m, bool project = false);
  void map_space_dims(const Dims& pfunc);
  void remove_space_dims(const Index_Set& vars);
  void remove_space_dim(Var var) {
    remove_space_dims(Index_Set(var.id()));
  }
  template <typename Iter>
  void remove_space_dims(Iter first, Iter last) {
    remove_space_dims(Index_Set(first, last));
  }
  void remove_higher_space_dims(dim_type new_dim);

  void concatenate_assign(const F_Poly& y);
  void expand_space_dim(Var var, dim_type m);
  void fold_space_dims(const Index_Set& vars, Var dest);

  // Semantically const, but may affect syntactic representation.
  void minimize() const;

  /* Input-output */
  bool ascii_load(std::istream& s);
  void ascii_dump(std::ostream& s) const;

  static F_Poly from_poly(const Poly& ph);
  Poly to_poly() const;

private:
  bool is_block_dim(dim_type i) const { return itvs[i].is_empty(); }
  bool is_itv_dim(dim_type i) const { return not is_block_dim(i); }
  bool are_itv_dims(const Block& b) const {
    return all_of(b, [this](dim_type i) { return is_itv_dim(i); });
  }

  bool is_boxable(const Factor& f) const {
    // TODO: if Itv will support NNC, remove check for closure.
    return f.space_dim() == 1 && f.is_topologically_closed();
  }

  // Moves i from itvs to blocks, adding it to existing block.
  void itv_to_block(dim_type i, dim_type block_idx);

  // Moves factor to itv componet, assuming it is "boxable".
  void block_to_itv(dim_type block_idx);
  // Moves boxable factors to itv component.
  void blocks_to_itvs();

  Cons itvs_to_cons() const;

  // Note: after merge, if b.size() > 1, then this will contain
  // a block having same dims as b, but possibly ordered differently.
  dim_type merge(const Block& b);

  void erase_empty_block(dim_type block_index);

  // Semantically const, but may affect syntactic representation.
  void sync(const Blocks& bs) const;
  void factorize(bool normalize = false) const;

  void normalize() const;

  void add_line_or_ray(const Gen& g);

}; // class F_Poly

NOTHROW_MOVES(F_Poly_Impl);
NOTHROW_MOVES(F_Poly);

namespace IO_Operators {

inline std::ostream&
operator<<(std::ostream& os, const F_Poly& ph) {
  ph.print(os);
  return os;
}

} // namespace IO_Operators

inline void
swap(F_Poly& x, F_Poly& y) noexcept { x.m_swap(y); }

inline bool
operator==(const F_Poly& x, const F_Poly& y) { return x.equals(y); }
inline bool
operator!=(const F_Poly& x, const F_Poly& y) { return !(x == y); }


namespace detail {

inline F_Poly_Impl::Block
extract_block(const Linear_Expr& e) {
  F_Poly_Impl::Block var_block;
  for (dim_type i : dim_range(e))
    if (e.get(Var(i)) != 0)
      var_block.push_back(i);
  return var_block;
}

inline F_Poly_Impl::Block
extract_block(const Con& c) {
  return extract_block(c.linear_expr());
}

inline F_Poly_Impl::Block
extract_block(const Gen& g) {
  assert(g.is_ray() || g.is_line());
  return extract_block(g.linear_expr());
}

inline bool
are_disjoint(const F_Poly_Impl::Block& b1,
             const F_Poly_Impl::Block& b2) {
  return std::find_first_of(b1.begin(), b1.end(),
                            b2.begin(), b2.end()) == b1.end();
}

} // namespace detail

template <typename Iter>
void
F_Poly::add_cons(Iter first, Iter last, bool do_normalize) {
  if (first == last || empty)
    return;
  if (do_normalize)
    factorize(true);

  // EXPERIMENTAL: sort cons according to number of affected blocks.
  using Iters = typename std::list<Iter>;
  std::map<dim_type, Iters> cs;

  auto count_affected
    = [this](const Block& b) {
        auto res = 0;
        for (auto i : b) {
          if (is_itv_dim(i))
            ++res;
        }
        for (const auto& b_i : blocks) {
          if (not detail::are_disjoint(b_i, b))
            ++res;
        }
        return res;
      };

  for (auto it = first; it != last; ++it) {
    Block b = detail::extract_block(*it);
    if (num_rows(b) <= 1) {
      // Be eager, as it won't require block merging
      add_con(*it);
    } else {
      auto num = count_affected(b);
      if (num == 1)
        // Be eager (ditto)
        add_con(*it);
      else
        // Be lazy, as it requires itv/block merging
        cs[num].push_back(it);
    }
  }

  // Add constraints ordered by number of block intersections.
  for (const auto& cluster : cs) {
    for (auto it : cluster.second) {
      add_con(*it);
    }
  }
  assert(check_inv());
}

template <typename Iter>
void
F_Poly::add_gens(Iter first, Iter last) {
  if (first == last)
    return;
  if (std::any_of(first, last, std::mem_fn(&Gen::is_point))) {
    Poly ph(space_dim(), Spec_Elem::EMPTY, topology());
    ph.add_gens(first, last);
    F_Poly y = from_poly(ph);
    y.factorize();
    poly_hull_assign(y);
    return;
  }

  // No point in the generators.
  assert(!empty);
  // Eagerly add all lines and rays.
  bool no_closure_points = true;
  for (auto it = first; it != last; ++it) {
    if (it->is_closure_point())
      no_closure_points = false;
    else
      add_line_or_ray(*it);
  }
  if (no_closure_points)
    return;

  // There were also closure points (but no points).
  auto& x = *this;
  Poly x_ph = x.to_poly();
  for (auto it = first; it != last; ++it) {
    if (it->is_closure_point())
      x_ph.add_gen(*first);
  }
  x = from_poly(x_ph);
  x.factorize();
}

} // namespace pplite

#endif // !defined(pplite_F_Poly_hh)
