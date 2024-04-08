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

#ifndef pplite_Dyn_Poly_hh
#define pplite_Dyn_Poly_hh 1

#include "globals.hh"
#include "Var.hh"
#include "Linear_Expr.hh"
#include "Bits.hh"
#include "Con.hh"
#include "Gen.hh"
#include "Poly_Rel.hh"
#include "Poly.hh"
#include "Abs_Poly.hh"

#include <memory>

namespace pplite {
namespace dynamic {

// An Abs_Poly wrapper, restoring C++ value semantics.
struct Dyn_Poly {
private:
  std::unique_ptr<Abs_Poly> ptr;
  // Private helper.
  explicit Dyn_Poly(Abs_Poly* abs_poly_ptr) : ptr(abs_poly_ptr) {}

public:
  explicit Dyn_Poly(dim_type d = 0,
                    Spec_Elem s = Spec_Elem::UNIVERSE,
                    Topol t = get_default_topology(),
                    Abs_Poly::Kind k = default_poly_kind)
    : ptr(make_poly(k, d, s, t)) {}

  // Allow for all possible orderings of the first three arguments.
  Dyn_Poly(dim_type d, Topol t, Spec_Elem s = Spec_Elem::UNIVERSE,
           Abs_Poly::Kind k = default_poly_kind)
    : Dyn_Poly(d, s, t, k) {}
  Dyn_Poly(Spec_Elem s, dim_type d, Topol t = get_default_topology(),
           Abs_Poly::Kind k = default_poly_kind)
    : Dyn_Poly(d, s, t, k) {}
  Dyn_Poly(Spec_Elem s, Topol t, dim_type d,
           Abs_Poly::Kind k = default_poly_kind)
    : Dyn_Poly(d, s, t, k) {}
  Dyn_Poly(Topol t, dim_type d, Spec_Elem s = Spec_Elem::UNIVERSE,
           Abs_Poly::Kind k = default_poly_kind)
    : Dyn_Poly(d, s, t, k) {}
  Dyn_Poly(Topol t, Spec_Elem s, dim_type d,
           Abs_Poly::Kind k = default_poly_kind)
    : Dyn_Poly(d, s, t, k) {}

  // Copy ctor and copy assign: reimplement with deep copy.
  Dyn_Poly(const Dyn_Poly& y) : ptr(y.ptr->clone()) {}
  Dyn_Poly& operator=(const Dyn_Poly& y) {
    Dyn_Poly temp = y;
    m_swap(temp);
    return *this;
  }
  // Move ctor, move assign and dtor: default.
  Dyn_Poly(Dyn_Poly&& y) = default;
  Dyn_Poly& operator=(Dyn_Poly&& y) = default;
  ~Dyn_Poly() = default;

  void minimize() const { ptr->minimize(); }

  Abs_Poly::Kind poly_kind() const { return ptr->poly_kind(); }
  bool is_necessarily_closed() const { return ptr->is_necessarily_closed(); }
  bool is_disjunctive() const { return ptr->is_disjunctive(); }
  bool check_inv() const { return ptr != nullptr && ptr->check_inv(); }

  bool is_empty() const { return ptr->is_empty(); }
  bool is_universe() const { return ptr->is_universe(); }
  bool is_bounded() const { return ptr->is_bounded(); }
  bool is_bounded_expr(bool from_below, const Linear_Expr& expr) const {
    return ptr->is_bounded_expr(from_below, expr);
  }
  bool is_topologically_closed() const {
    return ptr->is_topologically_closed();
  }
  bool boxed_contains(const Dyn_Poly& y) const {
    return ptr->boxed_contains(*y.ptr);
  }
  bool constrains(Var v) const { return ptr->constrains(v); }
  bool contains(const Dyn_Poly& y) const { return ptr->contains(*y.ptr); }
  bool equals(const Dyn_Poly& y) const { return ptr->equals(*y.ptr); }
  bool is_disjoint_from(const Dyn_Poly& y) const {
    return ptr->is_disjoint_from(*y.ptr);
  }

  Topol topology() const { return ptr->topology(); }
  dim_type space_dim() const { return ptr->space_dim(); }
  dim_type affine_dim() const { return ptr->affine_dim(); }
  dim_type num_min_cons() const { return ptr->num_min_cons(); }
  dim_type num_min_gens() const { return ptr->num_min_gens(); }
  Cons copy_cons() const { return ptr->copy_cons(); }
  Gens copy_gens() const { return ptr->copy_gens(); }
  Poly_Con_Rel relation_with(const Con& c) const {
    return ptr->relation_with(c);
  }
  Poly_Gen_Rel relation_with(const Gen& g) const {
    return ptr->relation_with(g);
  }
  BBox get_bounding_box() const { return ptr->get_bounding_box(); }
  bool min(const Affine_Expr& ae, Rational& value,
           bool* included_ptr = nullptr, Gen* g_ptr = nullptr) const {
    return ptr->min(ae, value, included_ptr, g_ptr);
  }
  bool max(const Affine_Expr& ae, Rational& value,
           bool* included_ptr = nullptr, Gen* g_ptr = nullptr) const {
    return ptr->max(ae, value, included_ptr, g_ptr);
  }
  Itv get_bounds(Var var) const { return ptr->get_bounds(var); }
  Itv get_bounds(const Affine_Expr& ae) const { return ptr->get_bounds(ae); }
  Index_Set get_unconstrained() const { return ptr->get_unconstrained(); }
  size_t hash() const { return ptr->hash(); }
  void print(std::ostream& os) const { ptr->print(os); }
  size_t get_memory_in_bytes() const { return ptr->get_memory_in_bytes(); }

  using Cons_Proxy = Abs_Poly::Cons_Proxy;
  Cons_Proxy cons() const { return ptr->cons(); }
  Cons_Proxy normalized_cons() const { return ptr->normalized_cons(); }
  using Gens_Proxy = Abs_Poly::Gens_Proxy;
  Gens_Proxy gens() const { return ptr->gens(); }

  void collapse(dim_type n) { ptr->collapse(n); }
  dim_type num_disjuncts() const { return ptr->num_disjuncts(); }
  Cons_Proxy disjunct_cons(dim_type n) const { return ptr->disjunct_cons(n); }
  bool geom_covers(const Dyn_Poly& y) const { return ptr->geom_covers(*y.ptr); }

  // Modifiers
  void m_swap(Dyn_Poly& y) noexcept { std::swap(ptr, y.ptr); }
  void set_empty() { ptr->set_empty(); }
  void set_universe() { ptr->set_universe(); }
  void set_topology(Topol t) { ptr->set_topology(t); }
  void add_con(const Con& c) { ptr->add_con(c); }
  void add_cons(const Cons& cs) { ptr->add_cons(cs); }
  void add_gen(const Gen& g) { ptr->add_gen(g); }
  void add_gens(const Gens& gs) { ptr->add_gens(gs); }

  template <typename Iter>
  void add_cons(Iter first, Iter last) {
    ptr->add_cons(first, last);
  }

  void con_hull_assign(const Dyn_Poly& y, bool boxed = false) {
    ptr->con_hull_assign(*y.ptr, boxed);
  }
  void concatenate_assign(const Dyn_Poly& y) {
    ptr->concatenate_assign(*y.ptr);
  }
  void intersection_assign(const Dyn_Poly& y) {
    ptr->intersection_assign(*y.ptr);
  }
  void join_assign(const Dyn_Poly& y) { ptr->join_assign(*y.ptr); }
  void poly_hull_assign(const Dyn_Poly& y) {
    ptr->poly_hull_assign(*y.ptr);
  }
  void poly_difference_assign(const Dyn_Poly& y) {
    ptr->poly_difference_assign(*y.ptr);
  }
  void time_elapse_assign(const Dyn_Poly& y) {
    ptr->time_elapse_assign(*y.ptr);
  }
  void topological_closure_assign() {
    ptr->topological_closure_assign();
  }

  void widening_assign(const Dyn_Poly& y,
                       Widen_Impl w_impl = get_widen_impl(),
                       Widen_Spec w_spec = get_widen_spec()) {
    ptr->widening_assign(*y.ptr, w_impl, w_spec);
  }
  void widening_assign(const Dyn_Poly& y, const Cons& cs,
                       Widen_Impl w_impl = get_widen_impl(),
                       Widen_Spec w_spec = get_widen_spec()) {
    ptr->widening_assign(*y.ptr, cs, w_impl, w_spec);
  }
  Dyn_Poly split(const Con& con, Topol topol) {
    return Dyn_Poly(ptr->split(con, topol));
  }
  Dyn_Poly split(const Con& c) { return Dyn_Poly(ptr->split(c)); }
  Dyn_Poly integral_split(const Con& c) {
    return Dyn_Poly(ptr->integral_split(c));
  }

  void
  affine_image(Var var, const Linear_Expr& expr,
               const Integer& inhomo = Integer::zero(),
               const Integer& den = Integer::one()) {
    ptr->affine_image(var, expr, inhomo, den);
  }
  void
  affine_preimage(Var var, const Linear_Expr& expr,
                  const Integer& inhomo = Integer::zero(),
                  const Integer& den = Integer::one()) {
    ptr->affine_preimage(var, expr, inhomo, den);
  }
  void
  parallel_affine_image(const Vars& vars,
                        const Linear_Exprs& exprs,
                        const Integers& inhomos,
                        const Integers& dens) {
    ptr->parallel_affine_image(vars, exprs, inhomos, dens);
  }

  void unconstrain(Var var) { ptr->unconstrain(var); }
  void unconstrain(const Index_Set& vars) { ptr->unconstrain(vars); }
  void add_space_dims(dim_type d, bool project = false) {
    ptr->add_space_dims(d, project);
  }
  void map_space_dims(const Dims& pfunc) { ptr->map_space_dims(pfunc); }
  void remove_space_dim(Var var) { ptr->remove_space_dim(var); }
  void remove_space_dims(const Index_Set& vars) {
    ptr->remove_space_dims(vars);
  }
  template <typename Iter>
  void remove_space_dims(Iter first, Iter last) {
    ptr->remove_space_dims(first, last);
  }
  void remove_higher_space_dims(dim_type new_dim) {
    ptr->remove_higher_space_dims(new_dim);
  }
  void expand_space_dim(Var var, dim_type m) {
    ptr->expand_space_dim(var, m);
  }
  void fold_space_dims(const Index_Set& vars, Var dest) {
    ptr->fold_space_dims(vars, dest);
  }

  /* Input-output */
  void ascii_dump(std::ostream& s) const { ptr->ascii_dump(s); }

}; // Dyn_Poly

inline void
swap(Dyn_Poly& x, Dyn_Poly& y) noexcept {
  x.m_swap(y);
}

inline bool
operator==(const Dyn_Poly& x, const Dyn_Poly& y) { return x.equals(y); }
inline bool
operator!=(const Dyn_Poly& x, const Dyn_Poly& y) { return !(x == y); }

} // namespace dynamic

namespace IO_Operators {

inline std::ostream&
operator<<(std::ostream& os, const dynamic::Dyn_Poly& x) {
  x.print(os);
  return os;
}

} // namespace IO_Operators

} // namespace pplite

#endif // pplite_Dyn_Poly_hh
