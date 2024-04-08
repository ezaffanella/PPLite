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

#ifndef pplite_PolySet_hh
#define pplite_PolySet_hh 1

#include "globals.hh"
#include "utils.hh"
#include "Bits.hh"
#include "BBox.hh"
#include "Con.hh"
#include "Gen.hh"
#include "Integer.hh"
#include "Itv.hh"
#include "Linear_Expr.hh"
#include "Poly_Rel.hh"
#include "Rational.hh"
#include "Var.hh"
#include "Poly.hh"
#include "B_Poly.hh"
#include "F_Poly.hh"
// #include "Poly_min.hh"
// #include "Poly_templ.hh"

#include <cassert>
#include <functional>
#include <iostream>
#include <list>
#include <map>

namespace pplite {

template <typename PH>
class PolySet {
public:
  using Disj = PH;
  using Seq = std::list<Disj>;
  struct Impl {
    dim_type dim;
    Topol topol;
    Seq seq;
    mutable bool reduced;
  };

private:
  Impl impl_;

public:
  Impl& impl() { return impl_; }
  const Impl& impl() const { return impl_; }
  Seq& seq() { return impl_.seq; }
  const Seq& seq() const { return impl_.seq; }

  bool is_reduced() const { return impl_.reduced; }
  void set_reduced() const { impl_.reduced = true; }
  void clear_reduced() const { impl_.reduced = false; }

  explicit
  PolySet(dim_type d = 0,
          Spec_Elem s = Spec_Elem::UNIVERSE,
          Topol t = get_default_topology())
    : impl_{d, t, Seq(), true} {
    if (s == Spec_Elem::UNIVERSE)
      seq().emplace_back(d, s, t);
    assert(check_inv());
  }

  PolySet(dim_type d, Topol t, Spec_Elem s = Spec_Elem::UNIVERSE)
    : PolySet(d, s, t) {}
  PolySet(Spec_Elem s, dim_type d, Topol t = get_default_topology())
    : PolySet(d, s, t) {}
  PolySet(Spec_Elem s, Topol t, dim_type d)
    : PolySet(d, s, t) {}
  PolySet(Topol t, dim_type d, Spec_Elem s = Spec_Elem::UNIVERSE)
    : PolySet(d, s, t) {}
  PolySet(Topol t, Spec_Elem s, dim_type d)
    : PolySet(d, s, t) {}

  explicit PolySet(Disj ph)
    : PolySet(ph.space_dim(), ph.topology(), Spec_Elem::EMPTY) {
    add_disjunct(std::move(ph));
    set_reduced();
    assert(check_inv());
  }

  PolySet(const PolySet& y) = default;
  PolySet& operator=(const PolySet& y) = default;
  PolySet(PolySet&& y) = default;
  PolySet& operator=(PolySet&& y) = default;
  ~PolySet() = default;

  /* Predicates */

  bool check_inv() const;

  bool is_necessarily_closed() const {
    return impl().topol == Topol::CLOSED;
  }

  bool is_empty() const {
    empty_reduce();
    return seq().empty();
  }

  bool is_minimized() const {
    return is_reduced()
      && all_of(seq(), std::mem_fn(&Disj::is_minimized));
  }

  bool is_universe() const {
    if (is_reduced())
      return has_single_disjunct() && seq().front().is_universe();
    auto res = any_of(seq(), std::mem_fn(&Disj::is_universe));
    if (res && not has_single_disjunct())
      const_cast<PolySet&>(*this).set_universe();
    return res;
  }

  bool is_topologically_closed() const {
    if (is_necessarily_closed())
      return true;
    omega_reduce();
    return all_of(seq(), std::mem_fn(&Disj::is_topologically_closed));
  }

  bool is_bounded() const {
    return all_of(seq(), std::mem_fn(&Disj::is_bounded));
  }

  bool is_bounded_expr(bool from_below, const Linear_Expr& expr) const {
    return all_of(seq(),
                  [from_below, &expr] (const Disj& ph) {
                    return ph.is_bounded_expr(from_below, expr);
                  });
  }

  bool constrains(Var var) const {
    assert(var.space_dim() <= space_dim());
    omega_reduce();
    return is_empty()
      || any_of(seq(), [var](const Disj& d) { return d.constrains(var); });
  }

  bool equals(const PolySet& y) const;

  bool contains(const PolySet& y) const {
    y.omega_reduce();
    if (y.is_empty())
      return true;
    omega_reduce();
    return all_of(y.seq(), [this](const Disj& d) { return contains(d); });
  }

  bool strictly_contains(const PolySet& y) const {
    return not equals(y) && contains(y);
  }

  bool is_disjoint_from(const PolySet& y) const {
    const auto& x = *this;
    x.omega_reduce();
    y.omega_reduce();
    if (x.is_empty() || y.is_empty())
      return true;
    return all_of(y.seq(),
                  [this](const Disj& d) {
                    return is_disjoint_from(d);
                  });
  }

  bool boxed_contains(const PolySet& y) const {
    // FIXME: unoptimized.
    return contains(y);
  }

  // Helpers.
  typename Seq::const_iterator begin() const { return seq().begin(); }
  typename Seq::const_iterator end() const { return seq().end(); }
  size_t size() const { return seq().size(); }
  bool has_single_disjunct() const {
    return begin() != end() && (++begin() == end());
  }
  bool contains(const Disj& dy) const {
    if (dy.is_empty())
      return true;
    return any_of(seq(), [&dy](const Disj& dx) { return dx.contains(dy); });
  }
  bool is_disjoint_from(const Disj& dy) const {
    if (dy.is_empty())
      return true;
    return all_of(seq(),
                  [&dy](const Disj& dx) { return dx.is_disjoint_from(dy); });
  }

  bool definitely_entails(const PolySet& y) const {
    assert(space_dim() == y.space_dim());
    const auto& x_seq = seq();
    const auto& y_seq = y.seq();
    for (const auto& xd : x_seq) {
      if (none_of(y_seq,
                  [&xd](const Disj& yd) { return yd.contains(xd); }))
        return false;
    }
    return true;
  }

  bool geom_covers(const PolySet& y) const;

  bool geom_equals(const PolySet& y) const {
    const auto& x = *this;
    return x.geom_covers(y) && y.geom_covers(x);
  }

  /* Queries */

  Topol topology() const { return impl().topol; }
  dim_type space_dim() const { return impl().dim; }
  dim_type affine_dim() const;

  /* Dummy impl */
  dim_type num_min_cons() const {
    omega_reduce();
    if (is_empty())
      return 1;
    dim_type res = 0;
    for (const auto& d : seq())
      res += d.num_min_cons();
    return res;
  }
  /* Dummy impl */
  dim_type num_min_gens() const {
    omega_reduce();
    dim_type res = 0;
    for (const auto& d : seq())
      res += d.num_min_gens();
    return res;
  }

  // (Partially) Unimplemented stuff.

  static const Disj& empty_disjunct() {
    static const Disj empty_d(0, Spec_Elem::EMPTY);
    return empty_d;
  }

  Cons copy_cons() const {
    if (has_single_disjunct())
      return begin()->copy_cons();
    if (is_empty())
      return empty_disjunct().copy_cons();
    PPLITE_UNIMPL;
  }

  Cons copy_normalized_cons() const {
    if (has_single_disjunct())
      return begin()->copy_normalized_cons();
    if (is_empty())
      return empty_disjunct().copy_normalized_cons();
    PPLITE_UNIMPL;
  }

  Gens copy_gens() const {
    if (has_single_disjunct())
      return begin()->copy_gens();
    if (is_empty())
      return empty_disjunct().copy_gens();
    PPLITE_UNIMPL;
  }

  using Cons_Proxy = typename Disj::Cons_Proxy;
  Cons_Proxy cons() const {
    if (has_single_disjunct())
      return begin()->cons();
    if (is_empty())
      return empty_disjunct().cons();
    PPLITE_UNIMPL;
  }

  Cons_Proxy normalized_cons() const {
    if (has_single_disjunct())
      return begin()->normalized_cons();
    if (is_empty())
      return empty_disjunct().normalized_cons();
    PPLITE_UNIMPL;
  }

  using Gens_Proxy = typename Disj::Gens_Proxy;
  Gens_Proxy gens() const {
    if (has_single_disjunct())
      return begin()->gens();
    if (is_empty())
      return empty_disjunct().gens();
    PPLITE_UNIMPL;
  }

  // Note: we do not force omega reduction here, results may vary.
  dim_type num_disjuncts() const { return size(); }
  Cons_Proxy disjunct_cons(dim_type n) const {
    assert(0 <= n && n < num_disjuncts());
    auto it = begin();
    std::advance(it, n);
    return it->cons();
  }

  Poly_Con_Rel relation_with(const Con& c) const;
  Poly_Gen_Rel relation_with(const Gen& g) const;

  bool min(const Affine_Expr& ae, Rational& value,
           bool* included_ptr = nullptr, Gen* g_ptr = nullptr) const;
  bool max(const Affine_Expr& ae, Rational& value,
           bool* included_ptr = nullptr, Gen* g_ptr = nullptr) const {
    bool res = min(-ae, value, included_ptr, g_ptr);
    if (res) neg_assign(value);
    return res;
  }

  Itv get_bounds(Var var) const {
    Itv res(Spec_Elem::EMPTY);
    for (const auto& d : seq()) {
      res.lub_assign(d.get_bounds(var));
      if (res.is_universe()) break;
    }
    return res;
  }

  Itv get_bounds(const Affine_Expr& ae) const {
    Itv res(Spec_Elem::EMPTY);
    for (const auto& d : seq()) {
      res.lub_assign(d.get_bounds(ae));
      if (res.is_universe()) break;
    }
    return res;
  }

  Itv get_bounds(const Itv_Expr& ie) const {
    Itv res(Spec_Elem::EMPTY);
    for (const auto& d : seq()) {
      res.lub_assign(d.get_bounds(ie));
      if (res.is_universe()) break;
    }
    return res;
  }

  Index_Set get_unconstrained() const {
    Index_Set res;
    if (is_empty())
      return res;
    res.set_until(space_dim());
    for (const auto& d : seq()) {
      res &= d.get_unconstrained();
      if (res.empty())
        break;
    }
    return res;
  }

  BBox get_bounding_box() const {
    BBox res(space_dim(), Spec_Elem::EMPTY);
    for (const auto& d : seq())
      res.lub_assign(d.get_bounding_box());
    return res;
  }

  size_t hash() const { PPLITE_UNIMPL; }

  /* Modifiers */

  void m_swap(PolySet& y) noexcept {
    using std::swap;
    swap(impl().dim, y.impl().dim);
    swap(impl().topol, y.impl().topol);
    swap(impl().seq, y.impl().seq);
    swap(impl().reduced, y.impl().reduced);
  }

  void set_empty() {
    seq().clear();
    set_reduced();
    assert(check_inv());
  }
  void set_universe() {
    set_empty();
    seq().emplace_back(space_dim(), topology(), Spec_Elem::UNIVERSE);
    assert(check_inv());
  }
  void set_topology(Topol t);

  void add_con(const Con& c) { add_cons(&c, &c + 1); }
  void add_cons(const Cons& cs) { add_cons(cs.begin(), cs.end()); }
  template <typename Iter>
  void add_cons(Iter first, Iter last);

  void add_gen(const Gen& g) { add_gens(&g, &g + 1); }
  void add_gens(const Gens& gs) { add_gens(gs.begin(), gs.end()); }
  template <typename Iter>
  void add_gens(Iter first, Iter last);

  void topological_closure_assign();

  template <typename Iter>
  void unconstrain(Iter first, Iter last);
  void unconstrain(Var var) {
    auto d = var.id();
    unconstrain(&d, &d + 1);
  }
  void unconstrain(const Index_Set& vars) {
    unconstrain(vars.begin(), vars.end());
  }

  void intersection_assign(const PolySet& y);

  void affine_image(Var var, const Linear_Expr& expr,
                    const Integer& inhomo = Integer::zero(),
                    const Integer& den = Integer::one()) {
    for (auto& d : seq()) {
      d.affine_image(var, expr, inhomo, den);
      clear_reduced();
    }
    assert(check_inv());
  }
  void affine_preimage(Var var, const Linear_Expr& expr,
                       const Integer& inhomo = Integer::zero(),
                       const Integer& den = Integer::one()) {
    for (auto& d : seq()) {
      d.affine_preimage(var, expr, inhomo, den);
      clear_reduced();
    }
    assert(check_inv());
  }
  void
  parallel_affine_image(const Vars& vars,
                        const Linear_Exprs& exprs,
                        const Integers& inhomos, const Integers& dens) {
    for (auto& d : seq()) {
      d.parallel_affine_image(vars, exprs, inhomos, dens);
      clear_reduced();
    }
    assert(check_inv());
  }

  void widening_assign(const PolySet& y, const Cons* upto_ptr,
                       Widen_Impl w_impl, Widen_Spec w_spec);
  void widening_assign(const PolySet& y,
                       Widen_Impl w_impl = get_widen_impl(),
                       Widen_Spec w_spec = get_widen_spec()) {
    widening_assign(y, nullptr, w_impl, w_spec);
  }
  void widening_assign(const PolySet& y, const Cons& upto_cons,
                       Widen_Impl w_impl = get_widen_impl(),
                       Widen_Spec w_spec = get_widen_spec()) {
    widening_assign(y, &upto_cons, w_impl, w_spec);
  }

  void time_elapse_assign(const PolySet& y) {
    pairwise_apply_assign(y, std::mem_fn(&Disj::time_elapse_assign));
  }

private:
  PolySet split_aux_ineq(const Con& c, Topol t, bool integral);
  PolySet integral_split_aux_eq(const Con& c);
public:
  PolySet split(const Con& c, Topol t) { return split_aux_ineq(c, t, false); }
  PolySet split(const Con& c) { return split(c, topology()); }
  PolySet integral_split(const Con& c) {
    if (c.is_equality())
      return integral_split_aux_eq(c);
    else
      return split_aux_ineq(c, topology(), true);
  }

  /* Change of space dim */
  void add_space_dims(dim_type m, bool project = false);
  void map_space_dims(const Dims& pfunc);
  void remove_higher_space_dims(dim_type new_dim);
  template <typename Iter>
  void remove_space_dims(Iter first, Iter last);
  void remove_space_dims(const Index_Set& vars) {
    remove_space_dims(vars.begin(), vars.end());
  }
  void remove_space_dim(Var var) {
    auto dim = var.id();
    remove_space_dims(&dim, &dim + 1);
  }
  void concatenate_assign(const PolySet& y);
  void expand_space_dim(Var var, dim_type m);
  void fold_space_dims(const Index_Set& vars, Var dest);

  // Semantically const, but may affect syntactic representation.
  void minimize() const {
    omega_reduce();
    for (const auto& d : seq())
      d.minimize();
  }

  // These are specific of powersets.

  static void omega_reduce_seq(Seq& seq) {
    auto i = seq.begin();
    while (i != seq.end()) {
      if (i->is_empty())
        i = seq.erase(i);
      else {
        auto j = seq.begin();
        while (j != seq.end()) {
          if (i != j && i->contains(*j))
            j = seq.erase(j);
          else
            ++j;
        }
        ++i;
      }
    }
  }

  void empty_reduce() const {
    if (is_reduced())
      return;
    auto& x_seq = const_cast<Seq&>(seq());
    x_seq.remove_if(std::mem_fn(&Disj::is_empty));
  }

  void omega_reduce() const {
    if (is_reduced())
      return;
    auto& x_seq = const_cast<Seq&>(seq());
    omega_reduce_seq(x_seq);
    set_reduced();
    assert(check_inv());
  }

  void add_disjunct(Disj d) {
    if (d.is_empty())
      return;
    seq().push_back(std::move(d));
    clear_reduced();
  }

  template <typename BinOpAssign>
  void
  pairwise_apply_assign(const PolySet& y, BinOpAssign op) {
    auto& x = *this;
    x.omega_reduce();
    y.omega_reduce();
    if (y.has_single_disjunct()) {
      // Special case (avoid copies of x disjuncts)
      const auto& dy = y.seq().front();
      for (auto& dx : x.seq())
        op(dx, dy);
    } else {
      // Need copies of x disjuncts
      Seq new_seq;
      for (const auto& dy : y.seq()) {
        for (auto dx : x.seq()) {
          op(dx, dy);
          if (dx.is_empty())
            continue;
          new_seq.emplace_back(std::move(dx));
        }
      }
      x.seq() = std::move(new_seq);
    }
    x.clear_reduced();
    assert(x.check_inv());
  }

  void poly_difference_assign(const PolySet& y) {
    // FIXME: should we require y to have a single disjunct?
    // Note: this is NOT difference_assign, so we collapse
    // into a single polyhedron (questionable!).
    collapse(1);
    if (is_empty() || y.is_empty())
      return;
    auto& dx = *seq().begin();
    auto it = y.begin();
    auto dy = *it;
    for (++it; it != y.end(); ++it)
      dy.poly_hull_assign(*it);
    dx.poly_difference_assign(dy);
    if (dx.is_empty())
      set_empty();
    assert(check_inv());
  }

  void poly_hull_assign(const PolySet& y) {
    collapse(1);
    if (is_empty()) {
      *this = y;
      collapse(1);
    } else {
      auto& dx = *seq().begin();
      for (const auto& dy : y)
        dx.poly_hull_assign(dy);
    }
    assert(check_inv());
  }

  void con_hull_assign(const PolySet& y, bool boxed = false) {
    auto& x = *this;
    x.omega_reduce();
    y.omega_reduce();
    SequenceAdapter<Disj> args;
    auto not_empty = [] (const Disj& d) { d.minimize(); return !d.is_empty(); };
    args.append(x.begin(), x.end(), not_empty);
    args.append(y.begin(), y.end(), not_empty);
    auto d = con_hull(args.begin(), args.end(), boxed);
    if (d.is_empty())
      x.set_empty();
    else
      x = PolySet(std::move(d));
    assert(check_inv());
  }

  void meet_assign(const PolySet& y) { intersection_assign(y); }
  void join_assign(const PolySet& y);
  void difference_assign(const PolySet& y);
  void collapse(dim_type n);

  /* Input-output */
  void print(std::ostream& s) const {
    if (is_empty()) {
      s << "empty";
      return;
    }
    auto it = begin();
    s << "{ ";
    it->print(s);
    s << " }";
    auto last = end();
    for (++it; it != last; ++it) {
      s << " || { ";
      it->print(s);
      s << " }";
    }
  }
  bool ascii_load(std::istream& s);
  void ascii_dump(std::ostream& s) const;

private:
  using DS_Pair = std::pair<Disj, Seq>;
  static DS_Pair linear_partition(const Disj& p, const Disj& q);
  static void linear_partition_aux(const Con& c, DS_Pair& p);
  static bool check_containment(const Disj& d, const Seq& seq);
  static bool check_omega_reduced(const Seq& seq);

}; // class PolySet

// Explicit template instantiation (declaration).
extern template class PolySet<B_Poly>;
extern template class PolySet<BF_Poly>;

using P_Set = PolySet<B_Poly>;
using FP_Set = PolySet<BF_Poly>;

NOTHROW_MOVES(P_Set);
NOTHROW_MOVES(FP_Set);

namespace IO_Operators {

template <typename PH>
std::ostream& operator<<(std::ostream& s, const PolySet<PH>& x);

} // namespace IO_Operators

template <typename PH>
inline void
swap(PolySet<PH>& x, PolySet<PH>& y) noexcept { x.m_swap(y); }

template <typename PH>
inline bool
operator==(const PolySet<PH>& x, const PolySet<PH>& y) { return x.equals(y); }
template <typename PH>
inline bool
operator!=(const PolySet<PH>& x, const PolySet<PH>& y) { return !(x == y); }

} // namespace pplite

#include "PolySet_templ.hh"

#endif // !defined(pplite_PolySet_hh)
