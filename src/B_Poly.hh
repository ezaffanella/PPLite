/* PPLite: a lightweight library for convex polyhedra derived from PPL.
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

#ifndef pplite_B_Poly_hh
#define pplite_B_Poly_hh 1

#include "globals.hh"
#include "BBox.hh"
#include "Poly.hh"
#include "Poly_min.hh"
#include "Poly_templ.hh"
#include "F_Poly.hh"

#include <cassert>
#include <iostream>
#include <iterator>
#include <numeric>
#include <vector>

namespace pplite {

template <typename PH>
class B_Wrap {
public:
  using poly_type = PH;
  using bbox_type = Box<true>; // Note: this keeps pseudo-volume info.
  using bbox_ptr = std::unique_ptr<bbox_type>;

  bool has_valid_bbox() const { return bbp != nullptr; }
  const poly_type& impl_poly() const { return ph; }
  const bbox_type& impl_bbox() const {
    assert(has_valid_bbox());
    return *bbp;
  }

private:
  poly_type ph;
  mutable bbox_ptr bbp;

  void reset_bbox(bbox_type* ptr = nullptr) const { bbp.reset(ptr); }
  void ensure_bbox() const {
    if (has_valid_bbox())
      assert(bbp->equals(bbox_type(ph.get_bounding_box())));
    else
      reset_bbox(new bbox_type(ph.get_bounding_box()));
  }

  static dim_type get_single_nonzero(const Linear_Expr& expr) {
    auto res = expr.last_nonzero();
    if (res == not_a_dim() || not expr.all_zeroes(0, res))
      return not_a_dim();
    return res;
  }
  static bool is_simple_expr(const Linear_Expr& expr) {
    auto res = expr.last_nonzero();
    if (res == not_a_dim() || not expr.all_zeroes(0, res))
      return not_a_dim();
    return res;
  }

public:
  explicit
  B_Wrap(dim_type d = 0, Spec_Elem s = Spec_Elem::UNIVERSE,
         Topol t = get_default_topology())
    : ph(d, s, t), bbp(new bbox_type(d, s)) {
    assert(check_inv());
  }
  B_Wrap(dim_type d, Topol t, Spec_Elem s = Spec_Elem::UNIVERSE)
    : B_Wrap(d, s, t) {}
  B_Wrap(Spec_Elem s, dim_type d, Topol t = get_default_topology())
    : B_Wrap(d, s, t) {}
  B_Wrap(Spec_Elem s, Topol t, dim_type d)
    : B_Wrap(d, s, t) {}
  B_Wrap(Topol t, dim_type d, Spec_Elem s = Spec_Elem::UNIVERSE)
    : B_Wrap(d, s, t) {}
  B_Wrap(Topol t, Spec_Elem s, dim_type d)
    : B_Wrap(d, s, t) {}

  explicit B_Wrap(poly_type&& ph)
    : ph(std::move(ph)) {}

  B_Wrap(const B_Wrap& y) : ph(y.ph) {
    if (y.bbp != nullptr)
      bbp.reset(new bbox_type(*y.bbp));
  }

  B_Wrap& operator=(const B_Wrap& y) {
    if (this != &y) {
      auto tmp = y;
      operator=(std::move(tmp));
    }
    return *this;
  }

  B_Wrap(B_Wrap&& y) noexcept = default;
  B_Wrap& operator=(B_Wrap&& y) noexcept = default;
  ~B_Wrap() = default;

  /* Predicates */
  bool check_inv() const {
    auto maybe_print = [](const char* s) {
#ifdef NDEBUG // Not in debugging mode: keep silent.
      (void) s;
#else // In debugging mode: be noisy.
      std::cerr << "B_Poly: " << s << "." << std::endl;
#endif
    };

    if (has_valid_bbox()) {
      if (not bbp->check_inv()) {
        maybe_print("bounding box is broken");
        return false;
      }
      if (ph.space_dim() != bbp->space_dim()) {
        maybe_print("different space dim");
        return false;
      }
      auto expected = bbox_type(ph.get_bounding_box());
      if (not bbp->equals(expected)) {
        maybe_print("bbox and poly are out of sync");
        return false;
      }
    }
    return ph.check_inv();
  }

  bool is_empty() const {
    return has_valid_bbox() ? bbp->is_empty() : ph.is_empty();
  }

  bool is_universe() const {
    if (has_valid_bbox() && not bbp->is_universe())
      return false;
    return ph.is_universe();
  }

  bool is_minimized() const { return ph.is_minimized(); }
  bool is_necessarily_closed() const { return ph.is_necessarily_closed(); }
  bool is_topologically_closed() const { return ph.is_topologically_closed(); }

  bool is_bounded() const {
    return has_valid_bbox() ? bbp->is_bounded() : ph.is_bounded();
  }
  bool is_bounded_expr(bool from_below, const Linear_Expr& expr) const {
    if (has_valid_bbox()) {
      if (bbp->is_bounded_expr(from_below, expr))
        return true;
      // false result is reliable for interval expressions
      if (is_interval_expr(expr))
        return false;
    }
    return ph.is_bounded_expr(from_below, expr);
  }

  bool constrains(Var var) const {
    if (has_valid_bbox() && bbp->constrains(var))
      return true;
    return ph.constrains(var);
  }

  bool equals(const B_Wrap& y) const {
    const auto& x = *this;
    assert(x.check_inv() && y.check_inv());
    x.ensure_bbox();
    y.ensure_bbox();
    return x.bbp->equals(*y.bbp)
      && x.ph.boxed_contains(y.ph)
      && y.ph.boxed_contains(x.ph);
  }

  bool contains(const B_Wrap& y) const {
    const auto& x = *this;
    x.ensure_bbox();
    y.ensure_bbox();
    return x.bbp->contains(*y.bbp) && x.ph.boxed_contains(y.ph);
  }

  bool strictly_contains(const B_Wrap& y) const {
    const auto& x = *this;
    x.ensure_bbox();
    y.ensure_bbox();
    return x.bbp->contains(*y.bbp) && x.ph.strictly_contains(y.ph);
  }

  bool is_disjoint_from(const B_Wrap& y) const {
    const auto& x = *this;
    x.ensure_bbox();
    y.ensure_bbox();
    return x.bbp->is_disjoint_from(*y.bbp)
      || x.ph.is_disjoint_from(y.ph);
  }

  bool boxed_contains(const B_Wrap& y) const {
    return ph.boxed_contains(y.ph);
  }

  /* Queries */
  Topol topology() const { return ph.topology(); }
  dim_type space_dim() const { return ph.space_dim(); }
  dim_type affine_dim() const { return ph.affine_dim(); }

  Poly_Con_Rel relation_with(const Con& c) const {
    return ph.relation_with(c);
  }
  Poly_Gen_Rel relation_with(const Gen& g) const {
    return ph.relation_with(g);
  }

  bool min(const Affine_Expr& ae, Rational& value,
           bool* included_ptr = nullptr, Gen* g_ptr = nullptr) const {
    // FIXME: improvable
    return ph.min(ae, value, included_ptr, g_ptr);
  }

  bool max(const Affine_Expr& ae, Rational& value,
           bool* included_ptr = nullptr, Gen* g_ptr = nullptr) const {
    if (not min(-ae, value, included_ptr, g_ptr))
      return false;
    neg_assign(value);
    return true;
  }

  BBox get_bounding_box() const {
    ensure_bbox();
    return BBox(*bbp);
  }

  Itv get_bounds(Var var) const {
    ensure_bbox();
    return bbp->get_bounds(var);
  }

  Itv get_bounds(const Affine_Expr& ae) const {
    // FIXME: improvable
    return ph.get_bounds(ae);
  }
  Itv get_bounds(const Itv_Expr& ie) const {
    // FIXME: improvable
    return ph.get_bounds(ie);
  }

  Index_Set get_unconstrained() const { return ph.get_unconstrained(); }

  Cons copy_cons() const { return ph.copy_cons(); }
  Gens copy_gens() const { return ph.copy_gens(); }
  Cons copy_normalized_cons() const { return ph.copy_normalized_cons(); }
  Gens_Info gens_info() const { return ph.gens_info(); }
  dim_type num_min_cons() const { return ph.num_min_cons(); }
  dim_type num_min_gens() const { return ph.num_min_gens(); }

  using Cons_Proxy = typename poly_type::Cons_Proxy;
  Cons_Proxy cons() const { return ph.cons(); }
  Cons_Proxy normalized_cons() const { return ph.normalized_cons(); }

  using Gens_Proxy = typename poly_type::Gens_Proxy;
  Gens_Proxy gens() const { return ph.gens(); }

  void collapse(dim_type n) {
    // bbox preserved
    ph.collapse(n);
  }
  dim_type num_disjuncts() const { return ph.num_disjuncts(); }
  Cons_Proxy disjunct_cons(dim_type n) const { return ph.disjunct_cons(n); }
  bool geom_covers(const B_Wrap& y) const {
    const auto& x = *this;
    if (x.has_valid_bbox() && y.has_valid_bbox()) {
      if (not x.bbp->contains(*y.bbp))
        return false;
    }
    return ph.geom_covers(y.ph);
  }

  /* Modifiers */

  void m_swap(B_Wrap& y) noexcept {
    using std::swap;
    swap(ph, y.ph);
    swap(bbp, y.bbp);
  }

  void set_empty() {
    if (has_valid_bbox())
      bbp->set_empty();
    ph.set_empty();
    assert(check_inv());
  }
  void set_universe() {
    if (has_valid_bbox())
      bbp->set_universe();
    ph.set_universe();
    assert(check_inv());
  }
  void set_topology(Topol t) {
    // bbox preserved
    ph.set_topology(t);
  }

  void add_con(Con c) {
    if (c.is_tautological())
      return;
    reset_bbox();
    ph.add_con(std::move(c));
  }
  void add_cons(Cons cs) {
    reset_bbox();
    ph.add_cons(std::move(cs));
  }
  template <typename Iter>
  void add_cons(Iter first, Iter last) {
    reset_bbox();
    ph.add_cons(first, last);
  }

  void add_gen(Gen g) {
    reset_bbox();
    ph.add_gen(std::move(g));
  }
  void add_gens(Gens gs) {
    reset_bbox();
    ph.add_gens(std::move(gs));
  }
  template <typename Iter>
  void add_gens(Iter first, Iter last) {
    reset_bbox();
    ph.add_gens(first, last);
  }

  void topological_closure_assign() {
    // bbox preserved
    ph.topological_closure_assign();
  }

  template <typename Iter>
  void unconstrain(Iter first, Iter last) {
    if (first == last)
      return;
    if (is_empty())
      return;
    if (has_valid_bbox())
      bbp->unconstrain(first, last);
    ph.unconstrain(first, last);
  }
  void unconstrain(Var var) {
    auto d = var.id();
    unconstrain(&d, &d + 1);
  }
  void unconstrain(const Index_Set& vars) {
    unconstrain(vars.begin(), vars.end());
  }

  void intersection_assign(const B_Wrap& y) {
    auto& x = *this;
    if (x.has_valid_bbox() && y.has_valid_bbox()) {
      auto& x_bb = *x.bbp;
      const auto& y_bb = *y.bbp;
      if (x_bb.is_empty())
        return;
      if (y_bb.is_empty()) {
        x.set_empty();
        return;
      }
      x_bb.glb_assign(y_bb);
      if (x_bb.is_empty()) {
        x.set_empty();
        return;
      }
    }
    x.reset_bbox();
    x.ph.intersection_assign(y.ph);
    assert(check_inv());
  }

  void join_assign(const B_Wrap& y) {
    if (y.is_empty())
      return;
    reset_bbox();
    ph.join_assign(y.ph);
  }

  void poly_difference_assign(const B_Wrap& y) {
    if (y.is_empty())
      return;
    reset_bbox();
    ph.poly_difference_assign(y.ph);
  }

  void poly_hull_assign(const B_Wrap& y) {
    if (y.is_empty())
      return;
    reset_bbox();
    ph.poly_hull_assign(y.ph);
  }

  void con_hull_assign(const B_Wrap& y, bool boxed = false) {
    if (y.is_empty())
      return;
    reset_bbox();
    ph.con_hull_assign(y.ph, boxed);
  }

  void affine_image(Var var, const Linear_Expr& expr,
                    const Integer& inhomo = Integer::zero(),
                    const Integer& den = Integer::one()) {
    // FIXME: improvable
    reset_bbox();
    ph.affine_image(var, expr, inhomo, den);
  }

  void affine_preimage(Var var, const Linear_Expr& expr,
                       const Integer& inhomo = Integer::zero(),
                       const Integer& den = Integer::one()) {
    reset_bbox();
    ph.affine_preimage(var, expr, inhomo, den);
  }

  void
  parallel_affine_image(const Vars& vars,
                        const Linear_Exprs& exprs,
                        const Integers& inhomos, const Integers& dens) {
    reset_bbox();
    ph.parallel_affine_image(vars, exprs, inhomos, dens);
  }

  void widening_assign(const B_Wrap& y,
                       Widen_Impl w_impl = get_widen_impl(),
                       Widen_Spec w_spec = get_widen_spec()) {
    reset_bbox();
    ph.widening_assign(y.ph, w_impl, w_spec);
  }
  void widening_assign(const B_Wrap& y, const Cons& upto_cons,
                       Widen_Impl w_impl = get_widen_impl(),
                       Widen_Spec w_spec = get_widen_spec()) {
    reset_bbox();
    ph.widening_assign(y.ph, upto_cons, w_impl, w_spec);
  }

  void time_elapse_assign(const B_Wrap& y) {
    reset_bbox();
    ph.time_elapse_assign(y.ph);
  }

  /* Change of space dim */
  void add_space_dims(dim_type m, bool project = false) {
    if (has_valid_bbox())
      bbp->add_space_dims(m, project);
    ph.add_space_dims(m, project);
  }

  void map_space_dims(const Dims& pfunc) {
    if (has_valid_bbox())
      bbp->map_space_dims(pfunc);
    ph.map_space_dims(pfunc);
  }

  void remove_higher_space_dims(dim_type new_dim) {
    if (has_valid_bbox())
      bbp->remove_higher_space_dims(new_dim);
    ph.remove_higher_space_dims(new_dim);
  }

  template <typename Iter>
  void remove_space_dims(Iter first, Iter last) {
    if (has_valid_bbox())
      bbp->remove_space_dims(first, last);
    ph.remove_space_dims(first, last);
  }
  void remove_space_dims(const Index_Set& vars) {
    remove_space_dims(vars.begin(), vars.end());
  }
  void remove_space_dim(Var var) {
    auto dim = var.id();
    remove_space_dims(&dim, &dim + 1);
  }

  void concatenate_assign(const B_Wrap& y) {
    auto& x = *this;
    if (x.has_valid_bbox() && y.has_valid_bbox())
      x.bbp->concatenate_assign(*y.bbp);
    else
      x.reset_bbox();
    x.ph.concatenate_assign(y.ph);
  }

  void expand_space_dim(Var var, dim_type m) {
    if (has_valid_bbox())
      bbp->expand_space_dim(var, m);
    ph.expand_space_dim(var, m);
  }

  void fold_space_dims(const Index_Set& vars, Var dest) {
    reset_bbox();
    ph.fold_space_dims(vars, dest);
  }

private:
  B_Wrap split_aux(const Con& c, Topol t, bool integral) {
    B_Wrap res(space_dim(), Spec_Elem::EMPTY, topology());
    if (ph.is_empty())
      return res;
    // FIXME: improvable? (check if box has empty intersection with c)
    res.ph = integral ? ph.integral_split(c) : ph.split(c, t);
    if (not res.ph.is_empty()) {
      reset_bbox();
      res.reset_bbox();
    }
    assert(check_inv() && res.check_inv());
    return res;
  }
public:
  B_Wrap split(const Con& c, Topol t) { return split_aux(c, t, false); }
  B_Wrap split(const Con& c) { return split(c, topology()); }
  B_Wrap integral_split(const Con& c) {
    return split_aux(c, topology(), true);
  }

  void minimize() const {
    ph.minimize();
    // minimization does not affect bbox
    assert(check_inv());
  }

  size_t hash() const { return ph.hash(); }

  /* Input-output */
  void print(std::ostream& s) const { ph.print(s); }
  void print() const { ph.print(); }
  bool ascii_load(std::istream& s) {
    // FIXME: also load bbox
    reset_bbox();
    return ph.ascii_load(s);
  }
  void ascii_dump(std::ostream& s) const {
    ph.ascii_dump(s);
    if (bbp == nullptr)
      s << "has_bbox = false\n";
    else {
      s << "has_bbox = true\n";
      bbp->ascii_dump(s);
    }
  }

}; // class B_Wrap

namespace IO_Operators {

template <typename PH>
inline std::ostream&
operator<<(std::ostream& s, const B_Wrap<PH>& ph) {
  ph.print(s);
  return s;
}

} // namespace IO_Operators

template <typename PH>
inline void
swap(B_Wrap<PH>& x, B_Wrap<PH>& y) noexcept { x.m_swap(y); }

template <typename PH>
inline bool
operator==(const B_Wrap<PH>& x, const B_Wrap<PH>& y) { return x.equals(y); }
template <typename PH>
inline bool
operator!=(const B_Wrap<PH>& x, const B_Wrap<PH>& y) { return !(x == y); }

// Explicit template instantiation (declaration).
extern template class B_Wrap<Poly>;
extern template class B_Wrap<F_Poly>;

using B_Poly = B_Wrap<Poly>;
using BF_Poly = B_Wrap<F_Poly>;

NOTHROW_MOVES(B_Poly);
NOTHROW_MOVES(BF_Poly);

} // namespace pplite

#endif // !defined(pplite_B_Poly_hh)
