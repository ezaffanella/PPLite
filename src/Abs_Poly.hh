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

#ifndef pplite_Abs_Poly_hh
#define pplite_Abs_Poly_hh 1

#include "globals.hh"
#include "mater_iterator.hh"
#include "Bits.hh"
#include "Con.hh"
#include "Gen.hh"
#include "Linear_Expr.hh"
#include "Poly_Rel.hh"
#include "Poly.hh"
#include "Var.hh"

namespace pplite {
namespace dynamic {

// A template for a polymorphic sequence
template <typename Value>
struct Abs_Sequence {
  using value_type = Value;

  Abs_Sequence() = default;
  virtual ~Abs_Sequence() {}

  Abs_Sequence(const Abs_Sequence&) = delete;
  Abs_Sequence(Abs_Sequence&&) = delete;
  Abs_Sequence& operator=(const Abs_Sequence&) = delete;
  Abs_Sequence& operator=(Abs_Sequence&&) = delete;

  // Note: size() may be less than end_pos() due to skippable values.
  virtual dim_type end_pos() const = 0;
  virtual dim_type size() const = 0;
  virtual bool is_skippable(dim_type pos) const = 0;
  virtual const value_type* get_value_ptr(dim_type pos) const = 0;
}; // Abs_Sequence

// The abstract class for polyhedra
// (note: an implementation detail; users should rather use Dyn_Poly)
struct Abs_Poly {
  Abs_Poly() = default;
  virtual ~Abs_Poly() = default;

  Abs_Poly(const Abs_Poly&) = delete;
  Abs_Poly(Abs_Poly&&) = delete;
  Abs_Poly& operator=(const Abs_Poly&) = delete;
  Abs_Poly& operator=(Abs_Poly&&) = delete;

  virtual Abs_Poly* clone() const = 0;

  virtual void minimize() const = 0;

  // Note: enum forward declaration
  enum class Kind : int;
  virtual Kind poly_kind() const = 0;
  virtual bool is_necessarily_closed() const = 0;
  virtual bool is_disjunctive() const = 0;
  virtual bool check_inv() const = 0;

  virtual bool is_empty() const = 0;
  virtual bool is_universe() const = 0;
  virtual bool is_bounded() const = 0;
  virtual bool is_bounded_expr(bool from_below,
                               const Linear_Expr& expr) const = 0;
  virtual bool is_topologically_closed() const = 0;
  virtual bool boxed_contains(const Abs_Poly& y) const = 0;
  virtual bool constrains(Var v) const = 0;
  virtual bool contains(const Abs_Poly& y) const = 0;
  virtual bool equals(const Abs_Poly& y) const = 0;
  virtual bool is_disjoint_from(const Abs_Poly& y) const = 0;

  virtual Topol topology() const = 0;
  virtual dim_type space_dim() const = 0;
  virtual dim_type affine_dim() const = 0;
  virtual dim_type num_min_cons() const = 0;
  virtual dim_type num_min_gens() const = 0;
  virtual Cons copy_cons() const = 0;
  virtual Gens copy_gens() const = 0;
  virtual Poly_Con_Rel relation_with(const Con& c) const = 0;
  virtual Poly_Gen_Rel relation_with(const Gen& g) const = 0;
  virtual BBox get_bounding_box() const = 0;
  virtual bool min(const Affine_Expr& ae, Rational& value,
                   bool* included_ptr = nullptr,
                   Gen* g_ptr = nullptr) const = 0;
  virtual bool max(const Affine_Expr& ae, Rational& value,
                   bool* included_ptr = nullptr,
                   Gen* g_ptr = nullptr) const = 0;
  virtual Itv get_bounds(Var var) const = 0;
  virtual Itv get_bounds(const Affine_Expr& ae) const = 0;
  virtual Itv get_bounds(const Itv_Expr& ie) const = 0;
  virtual Index_Set get_unconstrained() const = 0;
  virtual size_t hash() const = 0;
  virtual size_t get_memory_in_bytes() const = 0;

  // Proxies for polymorphic sequences (trying to minimize copies)
  template <typename Value>
  struct Proxy {
    using value_type = Value;
    using seq_type = Abs_Sequence<Value>;
    using seq_handle = std::unique_ptr<const seq_type>;
    seq_handle seq;

    explicit Proxy(const seq_type* seq_ptr) {
      seq.reset(seq_ptr);
    }

    dim_type size() const { return seq->size(); }
    dim_type end_pos() const { return seq->end_pos(); }
    bool is_skippable(dim_type pos) const { return seq->is_skippable(pos); }
    const value_type* get_value_ptr(dim_type pos) const {
      return seq->get_value_ptr(pos);
    }

    using const_iterator = Proxy_Iter<seq_type>;
    const_iterator begin() const { return const_iterator(seq.get(), false); }
    const_iterator end() const { return const_iterator(seq.get(), true); }
    const_iterator cbegin() const { return begin(); }
    const_iterator cend() const { return end(); }

  }; // Proxy

  using Cons_Proxy = Proxy<Con>;
  virtual Cons_Proxy cons() const = 0;
  virtual Cons_Proxy normalized_cons() const = 0;
  using Gens_Proxy = Proxy<Gen>;
  virtual Gens_Proxy gens() const = 0;

  // Extension for powerset-like abstract domains.
  virtual void collapse(dim_type n) = 0;
  virtual dim_type num_disjuncts() const = 0;
  virtual Cons_Proxy disjunct_cons(dim_type n) const = 0;
  virtual bool geom_covers(const Abs_Poly& y) const = 0;

  // Modifiers
  virtual void m_swap(Abs_Poly& y) noexcept = 0;
  virtual void set_empty() = 0;
  virtual void set_universe() = 0;
  virtual void set_topology(Topol t) = 0;
  virtual void add_con(const Con& c) = 0;
  virtual void add_cons(const Cons& cs) = 0;
  virtual void add_gen(const Gen& g) = 0;
  virtual void add_gens(const Gens& gs) = 0;

  template <typename Iter>
  void add_cons(Iter first, Iter last) {
    for ( ; first != last; ++first)
      add_con(*first);
  }

  virtual void con_hull_assign(const Abs_Poly& y, bool boxed = false) = 0;
  virtual void concatenate_assign(const Abs_Poly& y) = 0;
  virtual void intersection_assign(const Abs_Poly& y) = 0;
  virtual void join_assign(const Abs_Poly& y) = 0;
  virtual void poly_hull_assign(const Abs_Poly& y) = 0;
  virtual void poly_difference_assign(const Abs_Poly& y) = 0;
  virtual void time_elapse_assign(const Abs_Poly& y) = 0;
  virtual void topological_closure_assign() = 0;

  virtual void
  widening_assign(const Abs_Poly& y,
                  Widen_Impl w_impl = get_widen_impl(),
                  Widen_Spec w_spec = get_widen_spec()) = 0;
  virtual void
  widening_assign(const Abs_Poly& y, const Cons& cs,
                  Widen_Impl w_impl = get_widen_impl(),
                  Widen_Spec w_spec = get_widen_spec()) = 0;

  virtual Abs_Poly* split(const Con& con, Topol topol) = 0;
  virtual Abs_Poly* split(const Con& c) = 0;
  virtual Abs_Poly* integral_split(const Con& c) = 0;

  virtual void
  affine_image(Var var, const Linear_Expr& expr,
               const Integer& inhomo = Integer::zero(),
               const Integer& den = Integer::one()) = 0;
  virtual void
  affine_preimage(Var var, const Linear_Expr& expr,
                  const Integer& inhomo = Integer::zero(),
                  const Integer& den = Integer::one()) = 0;
  virtual void
  parallel_affine_image(const Vars& vars,
                        const Linear_Exprs& exprs,
                        const Integers& inhomos,
                        const Integers& dens) = 0;

  virtual void unconstrain(Var var) = 0;
  virtual void unconstrain(const Index_Set& vars) = 0;
  virtual void add_space_dims(dim_type d, bool project = false) = 0;
  virtual void map_space_dims(const Dims& pfunc) = 0;
  virtual void remove_space_dim(Var var) = 0;
  virtual void remove_space_dims(const Index_Set& vars) = 0;
  template <typename Iter>
  void remove_space_dims(Iter first, Iter last) {
    Index_Set vars(first, last);
    remove_space_dims(vars);
  }
  virtual void remove_higher_space_dims(dim_type new_dim) = 0;
  virtual void expand_space_dim(Var var, dim_type m) = 0;
  virtual void fold_space_dims(const Index_Set& vars, Var dest) = 0;

  /* Input-output */
  virtual void print(std::ostream& os) const = 0;
  virtual void ascii_dump(std::ostream& s) const = 0;

}; // Abs_Poly

Abs_Poly*
make_poly(Abs_Poly::Kind kind, dim_type d, Spec_Elem spec, Topol t);

Abs_Poly*
make_poly(const char* kind_name, dim_type d, Spec_Elem spec, Topol t);

} // namespace dynamic
} // namespace pplite

#endif // pplite_Abs_Poly_hh
