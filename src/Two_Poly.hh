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

#ifndef pplite_Two_Poly_hh
#define pplite_Two_Poly_hh 1

#include "globals.hh"

// For instantiation.
#include "Poly.hh"
#include "Poly_min.hh"
#include "Poly_templ.hh"
#include "PolySet.hh"

#include <cassert>
#include <iostream>

namespace pplite {

template <typename PH1, typename PH2>
class Two_Poly {
public:
  using Poly1 = PH1;
  using Poly2 = PH2;

  // Meant to be public.
  Poly1 ph1;
  Poly2 ph2;

  explicit
  Two_Poly(dim_type d = 0, Spec_Elem s = Spec_Elem::UNIVERSE,
           Topol t = get_default_topology())
    : ph1(d, s, t), ph2(d, s, t) { check_equiv(); }
  Two_Poly(dim_type d, Topol t, Spec_Elem s = Spec_Elem::UNIVERSE)
    : Two_Poly(d, s, t) {}
  Two_Poly(Spec_Elem s, dim_type d, Topol t = get_default_topology())
    : Two_Poly(d, s, t) {}
  Two_Poly(Spec_Elem s, Topol t, dim_type d)
    : Two_Poly(d, s, t) {}
  Two_Poly(Topol t, dim_type d, Spec_Elem s = Spec_Elem::UNIVERSE)
    : Two_Poly(d, s, t) {}
  Two_Poly(Topol t, Spec_Elem s, dim_type d)
    : Two_Poly(d, s, t) {}

  Two_Poly(const Two_Poly& y) = default;
  Two_Poly& operator=(const Two_Poly& y) = default;
  Two_Poly(Two_Poly&& y) noexcept = default;
  Two_Poly& operator=(Two_Poly&& y) noexcept = default;
  ~Two_Poly() = default;

  void check_equiv() const {
    assert(check_inv());
    // FIXME: this won't work on powersets (scan disjuncts?)
    static_assert((not is_PolySet<Poly1>::value &&
                   not is_PolySet<Poly2>::value),
                  "Two_Poly cannot handle powersets of polyhedra");
    Poly2 ph1_2(ph1.space_dim(), ph1.topology());
    ph1_2.add_cons(ph1.copy_cons());
    if (ph1_2 == ph2)
      return;
    using namespace IO_Operators;
    std::cerr << "=== Two_Poly: mismatch start ===\n";
    std::cerr << "=== Poly ph1 ===\n";
    std::cerr << ph1 << std::endl;
    std::cerr << "=== Poly ph2 ===\n";
    std::cerr << ph2 << std::endl;
    std::cerr << "=== Two_Poly mismatch end ===\n";
    abort();
  }

  void check_prop(bool prop) const {
    if (!prop) {
      using namespace IO_Operators;
      std::cerr << "=== Two_Poly: property does not hold ===\n";
      std::cerr << "=== Poly ph1 ===\n";
      std::cerr << ph1 << std::endl;
      std::cerr << "=== Poly ph2 ===\n";
      std::cerr << ph2 << std::endl;
      abort();
    }
  }

  bool check_inv() const {
    return ph1.check_inv() && ph2.check_inv();
  }

  /* Predicates */
  bool is_necessarily_closed() const {
    auto res1 = ph1.is_necessarily_closed();
    check_prop(res1 == ph2.is_necessarily_closed());
    return res1;
  }
  bool is_empty() const {
    auto res1 = ph1.is_empty();
    check_prop(res1 == ph2.is_empty());
    return res1;
  }
  bool is_minimized() const {
    auto res1 = ph1.is_minimized();
    check_prop(res1 == ph2.is_minimized());
    return res1;
  }
  bool is_universe() const {
    auto res1 = ph1.is_universe();
    check_prop(res1 == ph2.is_universe());
    return res1;
  }
  bool is_topologically_closed() const {
    auto res1 = ph1.is_topologically_closed();
    check_prop(res1 == ph2.is_topologically_closed());
    return res1;
  }
  bool is_bounded() const {
    auto res1 = ph1.is_bounded();
    check_prop(res1 == ph2.is_bounded());
    return res1;
  }
  bool is_bounded_expr(bool from_below, const Linear_Expr& expr) const {
    auto res1 = ph1.is_bounded_expr(from_below, expr);
    check_prop(res1 == ph2.is_bounded_expr(from_below, expr));
    return res1;
  }
  bool constrains(Var var) const {
    auto res1 = ph1.constrains(var);
    check_prop(res1 == ph2.constrains(var));
    return res1;
  }
  bool equals(const Two_Poly& y) const {
    auto res1 = ph1.equals(y.ph1);
    check_prop(res1 == ph2.equals(y.ph2));
    return res1;
  }
  bool contains(const Two_Poly& y) const {
    auto res1 = ph1.contains(y.ph1);
    check_prop(res1 == ph2.contains(y.ph2));
    return res1;
  }
  bool strictly_contains(const Two_Poly& y) const {
    auto res1 = ph1.strictly_contains(y.ph1);
    check_prop(res1 == ph2.strictly_contains(y.ph2));
    return res1;
  }
  bool is_disjoint_from(const Two_Poly& y) const {
    auto res1 = ph1.is_disjoint_from(y.ph1);
    check_prop(res1 == ph2.is_disjoint_from(y.ph2));
    return res1;
  }
  bool boxed_contains(const Two_Poly& y) const {
    auto res1 = ph1.boxed_contains(y.ph1);
    check_prop(res1 == ph2.boxed_contains(y.ph2));
    return res1;
  }

  /* Queries */
  Topol topology() const {
    auto res1 = ph1.topology();
    check_prop(res1 == ph2.topology());
    return res1;
  }
  dim_type space_dim() const {
    auto res1 = ph1.space_dim();
    check_prop(res1 == ph2.space_dim());
    return res1;
  }
  dim_type affine_dim() const {
    auto res1 = ph1.affine_dim();
    check_prop(res1 == ph2.affine_dim());
    return res1;
  }
  Poly_Con_Rel relation_with(const Con& c) const {
    auto res1 = ph1.relation_with(c);
    check_prop(res1 == ph2.relation_with(c));
    return res1;
  }
  Poly_Gen_Rel relation_with(const Gen& g) const {
    auto res1 = ph1.relation_with(g);
    check_prop(res1 == ph2.relation_with(g));
    return res1;
  }
  bool min(const Affine_Expr& ae, Rational& value,
           bool* included_ptr = nullptr, Gen* g_ptr = nullptr) const {
    auto value2 = value;
    bool included2;
    bool* included_ptr2 = included_ptr ? &included2 : nullptr;
    if (included_ptr)
      included2 = *included_ptr;
    Gen g2;
    Gen* g_ptr2 = g_ptr ? &g2 : nullptr;
    auto res1 = ph1.min(ae, value, included_ptr, g_ptr);
    auto res2 = ph2.min(ae, value2, included_ptr2, g_ptr2);
    check_prop(res1 == res2 && value == value2);
    if (included_ptr)
      check_prop(*included_ptr == *included_ptr2);
    // NOTE: we cannot pretend that the two domains actually
    // compute the same generator (i.e., in general *g_ptr != *g_ptr2).
    return res1;
  }
  bool max(const Affine_Expr& ae, Rational& value,
           bool* included_ptr = nullptr, Gen* g_ptr = nullptr) const {
    auto value2 = value;
    bool included2;
    bool* included_ptr2 = included_ptr ? &included2 : nullptr;
    if (included_ptr)
      included2 = *included_ptr;
    Gen g2;
    Gen* g_ptr2 = g_ptr ? &g2 : nullptr;
    auto res1 = ph1.max(ae, value, included_ptr, g_ptr);
    auto res2 = ph2.max(ae, value2, included_ptr2, g_ptr2);
    check_prop(res1 == res2 && value == value2);
    if (included_ptr)
      check_prop(*included_ptr == *included_ptr2);
    // NOTE: we cannot pretend that the two domains actually
    // compute the same generator (i.e., in general *g_ptr != *g_ptr2).
    return res1;
  }

  Itv get_bounds(Var var) const {
    auto res1 = ph1.get_bounds(var);
    auto res2 = ph2.get_bounds(var);
    check_prop(res1 == res2);
    return res1;
  }
  Itv get_bounds(const Affine_Expr& ae) const {
    auto res1 = ph1.get_bounds(ae);
    auto res2 = ph2.get_bounds(ae);
    check_prop(res1 == res2);
    return res1;
  }
  Itv get_bounds(const Itv_Expr& ie) const {
    auto res1 = ph1.get_bounds(ie);
    auto res2 = ph2.get_bounds(ie);
    check_prop(res1 == res2);
    return res1;
  }

  Index_Set get_unconstrained() const {
    auto res1 = ph1.get_unconstrained();
    auto res2 = ph2.get_unconstrained();
    check_prop(res1 == res2);
    return res1;
  }

  using Cons_Proxy = Cont_Proxy<Cons>;
  Cons_Proxy cons() const { return Cons_Proxy(copy_cons()); }
  using Gens_Proxy = Cont_Proxy<Gens>;
  Gens_Proxy gens() const { return Gens_Proxy(copy_gens()); }

  Cons copy_cons() const {
    auto res1 = ph1.copy_cons();
    // Have to force minimization to check for equality.
    Poly1 ph1_copy = ph1;
    ph1_copy.minimize();
    Poly2 ph2_copy = ph2;
    ph2_copy.minimize();
    auto res1_copy = ph2.copy_cons();
    auto res2_copy = ph2.copy_cons();
    make_set(res1_copy);
    make_set(res2_copy);
    check_prop(res1_copy == res2_copy);
    return res1;
  }
  Gens copy_gens() const {
    auto res1 = ph1.copy_gens();
    // Have to force minimization to check for equality.
    Poly1 ph1_copy = ph1;
    ph1_copy.minimize();
    Poly2 ph2_copy = ph2;
    ph2_copy.minimize();
    auto res1_copy = ph2.copy_gens();
    auto res2_copy = ph2.copy_gens();
    make_set(res1_copy);
    make_set(res2_copy);
    check_prop(res1_copy == res2_copy);
    return res1;
  }

  Cons_Proxy normalized_cons() const {
    return Cons_Proxy(copy_normalized_cons());
  }
  Cons copy_normalized_cons() const {
    auto res1 = ph1.copy_normalized_cons();
    // Have to force minimization to check for equality.
    Poly1 ph1_copy = ph1;
    auto res1_copy = ph2.copy_normalized_cons();
    Poly2 ph2_copy = ph2;
    auto res2_copy = ph2.copy_normalized_cons();
    make_set(res1_copy);
    make_set(res2_copy);
    check_prop(res1_copy == res2_copy);
    return res1;
  }

  void collapse(dim_type n) {
    ph1.collapse(n);
    ph2.collapse(n);
    check_equiv();
  }
  dim_type num_disjuncts() const {
    auto res1 = ph1.num_disjuncts();
    check_prop(res1 == ph2.num_disjuncts());
    return res1;
  }
  Cons_Proxy disjunct_cons(dim_type n) const {
    return ph1.disjunct_cons(n);
  }
  bool geom_covers(const Two_Poly& y) const {
    auto res1 = ph1.geom_covers(y.ph1);
    check_prop(res1 == ph2.geom_covers(y.ph2));
    return res1;
  }

  dim_type num_min_cons() const {
    auto res1 = ph1.num_min_cons();
    check_prop(res1 == ph2.num_min_cons());
    return res1;
  }
  dim_type num_min_gens() const {
    auto res1 = ph1.num_min_gens();
    check_prop(res1 == ph2.num_min_gens());
    return res1;
  }

  BBox get_bounding_box() const {
    auto res1 = ph1.get_bounding_box();
    check_prop(res1 == ph2.get_bounding_box());
    return res1;
  }

  /* Modifiers */

  void m_swap(Two_Poly& y) noexcept {
    using std::swap;
    swap(ph1, y.ph1);
    swap(ph2, y.ph2);
  }

  void set_empty() {
    ph1.set_empty();
    ph2.set_empty();
    check_equiv();
  }
  void set_universe() {
    ph1.set_universe();
    ph2.set_universe();
    check_equiv();
  }
  void set_topology(Topol t) {
    ph1.set_topology(t);
    ph2.set_topology(t);
    check_equiv();
  }

  void add_con(Con c) { add_cons(&c, &c + 1); }
  // Note: no point in passing move_iterators.
  void add_cons(Cons cs) { add_cons(cs.begin(), cs.end()); }
  template <typename Iter>
  void add_cons(Iter first, Iter last) {
    ph1.add_cons(first, last);
    ph2.add_cons(first, last);
    check_equiv();
  }

  void add_gen(Gen g) { add_gens(&g, &g + 1); }
  // Note: no point in passing move_iterators.
  void add_gens(Gens gs) { add_gens(gs.begin(), gs.end()); }
  template <typename Iter>
  void add_gens(Iter first, Iter last) {
    ph1.add_gens(first, last);
    ph2.add_gens(first, last);
    check_equiv();
  }

  Two_Poly split(const Con& c, Topol t) {
    Two_Poly res;
    res.ph1 = ph1.split(c, t);
    res.ph2 = ph2.split(c, t);
    check_equiv();
    res.check_equiv();
    return res;
  }
  Two_Poly split(const Con& c) { return split(c, topology()); }
  Two_Poly integral_split(const Con& c) {
    Two_Poly res;
    res.ph1 = ph1.integral_split(c);
    res.ph2 = ph2.integral_split(c);
    check_equiv();
    res.check_equiv();
    return res;
  }

  void topological_closure_assign() {
    ph1.topological_closure_assign();
    ph2.topological_closure_assign();
    check_equiv();
  }

  template <typename Iter>
  void unconstrain(Iter first, Iter last) {
    ph1.unconstrain(first, last);
    ph2.unconstrain(first, last);
    check_equiv();
  }
  void unconstrain(Var var) {
    auto d = var.id();
    unconstrain(&d, &d + 1);
  }
  void unconstrain(const Index_Set& vars) {
    unconstrain(vars.begin(), vars.end());
  }

  void intersection_assign(const Two_Poly& y) {
    ph1.intersection_assign(y.ph1);
    ph2.intersection_assign(y.ph2);
    check_equiv();
  }
  void join_assign(const Two_Poly& y) {
    ph1.join_assign(y.ph1);
    ph2.join_assign(y.ph2);
    check_equiv();
  }
  void con_hull_assign(const Two_Poly& y, bool boxed = false) {
    ph1.con_hull_assign(y.ph1, boxed);
    ph2.con_hull_assign(y.ph2, boxed);
    check_equiv();
  }
  void poly_hull_assign(const Two_Poly& y) {
    ph1.poly_hull_assign(y.ph1);
    ph2.poly_hull_assign(y.ph2);
    check_equiv();
  }
  void poly_difference_assign(const Two_Poly& y) {
    ph1.poly_difference_assign(y.ph1);
    ph2.poly_difference_assign(y.ph2);
    check_equiv();
  }

  void affine_image(Var var, const Linear_Expr& expr,
                    const Integer& inhomo = Integer::zero(),
                    const Integer& den = Integer::one()) {
    ph1.affine_image(var, expr, inhomo, den);
    ph2.affine_image(var, expr, inhomo, den);
    check_equiv();
  }
  void affine_preimage(Var var, const Linear_Expr& expr,
                       const Integer& inhomo = Integer::zero(),
                       const Integer& den = Integer::one()) {
    ph1.affine_preimage(var, expr, inhomo, den);
    ph2.affine_preimage(var, expr, inhomo, den);
    check_equiv();
  }
  void
  parallel_affine_image(const Vars& vars, const Linear_Exprs& exprs,
                        const Integers& inhomos, const Integers& dens) {
    ph1.parallel_affine_image(vars, exprs, inhomos, dens);
    ph2.parallel_affine_image(vars, exprs, inhomos, dens);
    check_equiv();
  }

  void widening_assign(const Two_Poly& y,
                       Widen_Impl w_impl = get_widen_impl(),
                       Widen_Spec w_spec = get_widen_spec()) {
    ph1.widening_assign(y.ph1, w_impl, w_spec);
    ph2.widening_assign(y.ph2, w_impl, w_spec);
    check_equiv();
  }
  void widening_assign(const Two_Poly& y, const Cons& upto_cons,
                       Widen_Impl w_impl = get_widen_impl(),
                       Widen_Spec w_spec = get_widen_spec()) {
    ph1.widening_assign(y.ph1, upto_cons, w_impl, w_spec);
    ph2.widening_assign(y.ph2, upto_cons, w_impl, w_spec);
    check_equiv();
  }

  void time_elapse_assign(const Two_Poly& y) {
    ph1.time_elapse_assign(y.ph1);
    ph2.time_elapse_assign(y.ph2);
    check_equiv();
  }

  /* Change of space dim */
  void add_space_dims(dim_type m, bool project = false) {
    ph1.add_space_dims(m, project);
    ph2.add_space_dims(m, project);
    check_equiv();
  }
  void map_space_dims(const Dims& pfunc) {
    ph1.map_space_dims(pfunc);
    ph2.map_space_dims(pfunc);
    check_equiv();
  }
  void remove_higher_space_dims(dim_type new_dim) {
    ph1.remove_higher_space_dims(new_dim);
    ph2.remove_higher_space_dims(new_dim);
    check_equiv();
  }
  template <typename Iter>
  void remove_space_dims(Iter first, Iter last) {
    ph1.remove_space_dims(first, last);
    ph2.remove_space_dims(first, last);
    check_equiv();
  }
  void remove_space_dims(const Index_Set& vars) {
    remove_space_dims(vars.begin(), vars.end());
  }
  void remove_space_dim(Var var) {
    auto dim = var.id();
    remove_space_dims(&dim, &dim + 1);
  }
  void concatenate_assign(const Two_Poly& y) {
    ph1.concatenate_assign(y.ph1);
    ph2.concatenate_assign(y.ph2);
    check_equiv();
  }
  void expand_space_dim(Var var, dim_type m) {
    ph1.expand_space_dim(var, m);
    ph2.expand_space_dim(var, m);
    check_equiv();
  }
  void fold_space_dims(const Index_Set& vars, Var dest) {
    ph1.fold_space_dims(vars, dest);
    ph2.fold_space_dims(vars, dest);
    check_equiv();
  }

  // Semantically const, but may affect syntactic representation.
  void minimize() const {
    ph1.minimize();
    ph2.minimize();
    check_equiv();
  }

  size_t hash() const {
    auto res1 = ph1.hash();
    check_equiv();
    // FIXME: hash is not checkable.
    return res1;
  }

  /* Input-output (here we only work on ph1) */
  void print(std::ostream& s) const { ph1.print(s); }
  void print() const { print(std::cout); }
  bool ascii_load(std::istream& s) {
    auto res = ph1.ascii_load(s);
    if (res) ph2 = Poly2(ph1);
    return res;
  }
  void ascii_dump(std::ostream& s) const {
    ph1.ascii_dump(s);
  }

private:
  template <typename T>
  static void
  make_set(std::vector<T>& vec) {
    std::sort(vec.begin(), vec.end(),
              [](const T& t1, const T& t2) {
                return compare(t1, t2) < 0;
              });
    vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
  }

}; // class Two_Poly

namespace IO_Operators {

template <typename PH1, typename PH2>
std::ostream& operator<<(std::ostream& s, const Two_Poly<PH1, PH2>& ph) {
  s << ph.ph1;
  return s;
}

} // namespace IO_Operators

template <typename PH1, typename PH2>
inline void
swap(Two_Poly<PH1, PH2>& x, Two_Poly<PH1, PH2>& y) noexcept { x.m_swap(y); }

template <typename PH1, typename PH2>
inline bool
operator==(const Two_Poly<PH1, PH2>& x, const Two_Poly<PH1, PH2>& y) {
  return x.equals(y);
}
template <typename PH1, typename PH2>
inline bool
operator!=(const Two_Poly<PH1, PH2>& x, const Two_Poly<PH1, PH2>& y) {
  return !(x == y);
}

} // namespace pplite

#endif // !defined(pplite_Two_Poly_hh)
