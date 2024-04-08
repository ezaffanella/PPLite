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

#ifndef pplite_Poly_hh
#define pplite_Poly_hh 1

#include "globals.hh"
#include "ascii_dump_load.hh"
#include "mater_iterator.hh"
#include "support_utils.hh"
#include "BBox.hh"
#include "Bits.hh"
#include "Con.hh"
#include "Gen.hh"
#include "Integer.hh"
#include "Itv.hh"
#include "Linear_Expr.hh"
#include "Poly_Rel.hh"
#include "Sat.hh"
#include "Var.hh"

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <tuple>
#include <utility>
#include <vector>

namespace pplite {

/* The implementation class: no encapsulation (everything public). */
struct Poly_Impl {
  enum class Status {
    EMPTY,
    MINIMIZED,
    PENDING
  };

  template <typename SK_Rows>
  struct Sys {
    using rows_type = SK_Rows;
    using row_type = typename SK_Rows::value_type;
    SK_Rows sing_rows;
    SK_Rows sk_rows;
    NS_Rows ns_rows;
    void clear() {
      sing_rows.clear();
      sk_rows.clear();
      ns_rows.clear();
    }
    bool empty() const {
      return sing_rows.empty() && sk_rows.empty();
    }
    void ascii_dump(std::ostream& s) const {
      s << "sing_rows " << sing_rows.size() << "\n";
      ascii_dump_all(s, sing_rows);
      s << "sk_rows " << sk_rows.size() << "\n";
      ascii_dump_all(s, sk_rows);
      s << "ns_rows " << ns_rows.size() << "\n";
      ascii_dump_all(s, ns_rows);
    }
    bool operator==(const Sys& y) const {
      return sing_rows == y.sing_rows
        && ns_rows == y.ns_rows
        && sk_rows == y.sk_rows;
    }
    void add(const row_type& new_elem) {
      if (detail::is_singular(new_elem))
        sing_rows.push_back(new_elem);
      else
        sk_rows.push_back(new_elem);
    }
  };

  NOTHROW_DEFAULT_AND_MOVES(Sys<Cons>);
  NOTHROW_DEFAULT_AND_MOVES(Sys<Gens>);

  // Basic status
  Topol topol;
  dim_type dim;
  Status status;
  // The DD pair.
  Sys<Cons> cs;
  Sys<Gens> gs;
  // Saturation info
  Sat sat_c;
  Sat sat_g;
  // Pending stuff
  Sys<Cons> cs_pending;
  Sys<Gens> gs_pending;

  // These methods are meant to be not trivially accessible from Poly.
  // They are anyway accessible using method Poly::impl().
  bool marked_empty() const { return status == Status::EMPTY; }
  bool marked_min() const { return status == Status::MINIMIZED; }
  bool has_pending() const { return status == Status::PENDING; }
  bool has_cs_pending() const {
    return has_pending() && !cs_pending.empty();
  }
  bool has_gs_pending() const {
    return has_pending() && !gs_pending.empty();
  }
  bool has_valid_cons() const { return !has_gs_pending(); }
  bool has_valid_gens() const { return !has_cs_pending(); }

  void set_status(Status s) { status = s; }

  void clear_all();
  void reinit_with_gens(dim_type, Topol, const Sys<Gens>&) = delete;
  void reinit_with_gens(dim_type d, Topol t, Sys<Gens>&& gens);

  void minimize() const; // Only logically const.
  static PPLITE_TLS dim_type minimize_filter_threshold;
  static dim_type get_minimize_filter_threshold() {
    return minimize_filter_threshold;
  }
  static void set_minimize_filter_threshold(dim_type threshold) {
    minimize_filter_threshold = threshold;
  }

  void ensure_valid_cons() const {
    if (!has_valid_cons())
      minimize();
    assert(check_inv(false));
  }
  void ensure_valid_gens() const {
    if (!has_valid_gens())
      minimize();
    assert(check_inv(false));
  }

  dim_type num_equals() const {
    assert(marked_empty() || marked_min());
    return cs.sing_rows.size();
  }
  dim_type num_lines() const {
    assert(marked_empty() || marked_min());
    return gs.sing_rows.size();
  }

  Gens_Info gens_info() const;
  dim_type num_min_cons() const;
  dim_type num_min_gens() const;

  std::tuple<const Sys<Cons>*, const Sys<Cons>*>
  get_sys_ptrs(const Sys<Cons>*) const {
    return std::make_tuple(&cs, &cs_pending);
  }
  std::tuple<const Sys<Gens>*, const Sys<Gens>*>
  get_sys_ptrs(const Sys<Gens>*) const {
    return std::make_tuple(&gs, &gs_pending);
  }

  template <typename Sys, typename Out, typename Pred>
  static Out
  materialize(const Sys& sys, const Sys& pending, Out out,
              bool cons_and_closed, Pred pred);

  // These methods are meant to be accessible from Poly.
  Poly_Impl(const Poly_Impl& y) = default;
  Poly_Impl(Poly_Impl&& y) = default;
  Poly_Impl& operator=(const Poly_Impl& y) = default;
  Poly_Impl& operator=(Poly_Impl&& y) = default;
  ~Poly_Impl() = default;

  explicit
  Poly_Impl(dim_type d = 0,
            Spec_Elem s = Spec_Elem::UNIVERSE,
            Topol t = get_default_topology());

  /* Predicates */

  bool check_inv(bool do_equals = true) const;
  bool is_necessarily_closed() const { return topol == Topol::CLOSED; }

  bool is_empty() const;
  bool is_minimized() const { return !has_pending(); }
  bool is_universe() const;
  bool is_topologically_closed() const;
  bool is_bounded() const;
  bool is_bounded_expr(bool from_below, const Linear_Expr& expr) const;

  bool constrains(Var var) const;
  bool equals(const Poly_Impl& y) const;
  bool contains(const Poly_Impl& y) const;
  bool strictly_contains(const Poly_Impl& y) const {
    const auto& x = *this;
    return x.contains(y) && !y.contains(x);
  }
  bool is_disjoint_from(const Poly_Impl& y) const;

  BBox get_bounding_box() const;
  bool boxed_contains(const Poly_Impl& y) const;

  /* Queries */

  Topol topology() const { return topol; }
  dim_type space_dim() const { return dim; }
  dim_type affine_dim() const;
  Poly_Con_Rel relation_with(const Con& c) const;
  Poly_Gen_Rel relation_with(const Gen& g) const;

  using Cons_Proxy = Mater_Sys<Sys<Cons>, Poly_Impl>;
  Cons_Proxy cons() const {
    ensure_valid_cons();
    return Cons_Proxy(this);
  }
  using Gens_Proxy = Mater_Sys<Sys<Gens>, Poly_Impl>;
  Gens_Proxy gens() const {
    ensure_valid_gens();
    return Gens_Proxy(this);
  }

  Cons_Proxy normalized_cons() const;
  Cons copy_normalized_cons() const {
    auto proxy = normalized_cons();
    return Cons(proxy.begin(), proxy.end());
  }

  Cons copy_cons() const {
    auto proxy = cons();
    return Cons(proxy.begin(), proxy.end());
  }
  Gens copy_gens() const {
    auto proxy = gens();
    return Gens(proxy.begin(), proxy.end());
  }

  bool min(const Affine_Expr& ae, Rational& value,
           bool* included_ptr = nullptr, Gen* g_ptr = nullptr) const;
  bool max(const Affine_Expr& ae, Rational& value,
           bool* included_ptr = nullptr, Gen* g_ptr = nullptr) const {
    bool res = min(-ae, value, included_ptr, g_ptr);
    if (res) neg_assign(value);
    return res;
  }

  Itv get_bounds(Var var) const;
  Itv get_bounds(const Affine_Expr& ae) const;
  Itv get_bounds(const Itv_Expr& ie) const;
  Index_Set get_unconstrained() const;

  size_t hash() const;
  void print(std::ostream& os) const;
  void print() const { print(std::cout); }

  void collapse(dim_type) { /* nothing to do */ }
  dim_type num_disjuncts() const { return is_empty() ? 0 : 1; }
  Cons_Proxy disjunct_cons(dim_type n) const {
    (void) n;
    assert(n == 0);
    return cons();
  }
  bool geom_covers(const Poly_Impl& y) const { return contains(y); }

  /* Modifiers */

  void m_swap(Poly_Impl& y) noexcept {
    using std::swap;
    swap(topol, y.topol);
    swap(dim, y.dim);
    swap(status, y.status);
    swap(cs, y.cs);
    swap(gs, y.gs);
    swap(sat_c, y.sat_c);
    swap(sat_g, y.sat_g);
    swap(cs_pending, y.cs_pending);
    swap(gs_pending, y.gs_pending);
  }

  void set_empty();
  void set_universe();
  void set_topology(Topol t);

  template <typename Iter>
  void add_cons(Iter first, Iter last);
  void add_con(Con c);
  void add_cons(Cons cs) {
    add_cons(std::make_move_iterator(cs.begin()),
             std::make_move_iterator(cs.end()));
  }

  template <typename Iter>
  void add_gens(Iter first, Iter last);
  void add_gen(Gen g);
  void add_gens(Gens gs) {
    add_gens(std::make_move_iterator(gs.begin()),
             std::make_move_iterator(gs.end()));
  }

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

  void intersection_assign(const Poly_Impl& y);
  void join_assign(const Poly_Impl& y) { poly_hull_assign(y); }
  void poly_hull_assign(const Poly_Impl& y);

  void affine_image(Var var, const Linear_Expr& expr,
                    const Integer& inhomo, const Integer& den);
  void affine_preimage(Var var, const Linear_Expr& expr,
                       const Integer& inhomo, const Integer& den);
  void
  parallel_affine_image(const Vars& vars,
                        const Linear_Exprs& exprs,
                        const Integers& inhomos, const Integers& dens);

  void widening_assign(const Poly_Impl& y, const Cons* upto_ptr,
                       Widen_Impl w_impl, Widen_Spec w_spec);

  void time_elapse_assign(const Poly_Impl& y);

  /* Change of space dim */

  void add_space_dims(dim_type m, bool project = false);

  // If the percentage of space dims to be removed is below this value,
  // efforts will be done to incrementally maintain the DD representation.
  static PPLITE_TLS int remove_space_dims_percentage;

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

  void map_space_dims(const Dims& pfunc);
  void concatenate_assign(const Poly_Impl& y);
  void expand_space_dim(Var var, dim_type m);
  void fold_space_dims(const Index_Set& vars, Var dest);

  Poly_Impl split(const Con& c, Topol t);
  Poly_Impl integral_split(const Con& c);

  /* Input-output */
  bool ascii_load(std::istream& s);
  void ascii_dump(std::ostream& s) const;

}; // struct Poly_Impl

inline void swap(Poly_Impl& x, Poly_Impl& y) noexcept { x.m_swap(y); }


/* The interface class: mild encapsulation, easy to overcome using impl(). */
class Poly : private Poly_Impl {
public:
  using Impl = Poly_Impl;
  explicit
  Poly(dim_type d = 0,
       Spec_Elem s = Spec_Elem::UNIVERSE,
       Topol t = get_default_topology())
    : Impl(d, s, t) {}

  // Allow for all possible orderings of the three arguments.
  Poly(dim_type d, Topol t, Spec_Elem s = Spec_Elem::UNIVERSE)
    : Impl(d, s, t) {}
  Poly(Spec_Elem s, dim_type d, Topol t = get_default_topology())
    : Impl(d, s, t) {}
  Poly(Spec_Elem s, Topol t, dim_type d)
    : Impl(d, s, t) {}
  Poly(Topol t, dim_type d, Spec_Elem s = Spec_Elem::UNIVERSE)
    : Impl(d, s, t) {}
  Poly(Topol t, Spec_Elem s, dim_type d)
    : Impl(d, s, t) {}

  Poly(const Poly& y) = default;
  Poly(Poly&& y) = default;
  Poly& operator=(const Poly& y) = default;
  Poly& operator=(Poly&& y) = default;
  ~Poly() = default;

  /* Accessors to low level implementation */
  Impl& impl() { return *this; }
  const Impl& impl() const { return *this; }

  /* Types */
  using Impl::Cons_Proxy;
  using Impl::Gens_Proxy;

  /* Predicates */
  bool is_necessarily_closed() const { return impl().is_necessarily_closed(); }
  bool check_inv(bool do_equals = true) const {
    return impl().check_inv(do_equals);
  }
  bool is_empty() const { return impl().is_empty(); }
  bool is_universe() const { return impl().is_universe(); }
  bool is_minimized() const { return impl().is_minimized(); }
  bool is_topologically_closed() const {
    return impl().is_topologically_closed();
  }
  bool is_bounded() const { return impl().is_bounded(); }
  bool is_bounded_expr(bool from_below, const Linear_Expr& expr) const {
    return impl().is_bounded_expr(from_below, expr);
  }
  bool constrains(Var var) const { return impl().constrains(var); }
  bool equals(const Poly& y) const {
    return impl().equals(y.impl());
  }
  bool contains(const Poly& y) const {
    return impl().contains(y.impl());
  }
  bool strictly_contains(const Poly& y) const {
    return impl().strictly_contains(y.impl());
  }
  bool is_disjoint_from(const Poly& y) const {
    return impl().is_disjoint_from(y.impl());
  }
  BBox get_bounding_box() const {
    return impl().get_bounding_box();
  }
  bool boxed_contains(const Poly& y) const {
    return impl().boxed_contains(y.impl());
  }

  /* Queries */
  Topol topology() const { return impl().topology(); }
  dim_type space_dim() const { return impl().space_dim(); }
  dim_type affine_dim() const { return impl().affine_dim(); }
  Poly_Con_Rel relation_with(const Con& c) const {
    return impl().relation_with(c);
  }
  Poly_Gen_Rel relation_with(const Gen& g) const {
    return impl().relation_with(g);
  }
  bool min(const Affine_Expr& ae, Rational& value,
           bool* included_ptr = nullptr, Gen* g_ptr = nullptr) const {
    return impl().min(ae, value, included_ptr, g_ptr);
  }
  bool max(const Affine_Expr& ae, Rational& value,
           bool* included_ptr = nullptr, Gen* g_ptr = nullptr) const {
    return impl().max(ae, value, included_ptr, g_ptr);
  }

  Itv get_bounds(Var var) const { return impl().get_bounds(var); }
  Itv get_bounds(const Affine_Expr& ae) const { return impl().get_bounds(ae); }
  Itv get_bounds(const Itv_Expr& ie) const { return impl().get_bounds(ie); }

  Index_Set get_unconstrained() const { return impl().get_unconstrained(); }

  size_t hash() const { return impl().hash(); }
  void print(std::ostream& os) const { impl().print(os); }
  void print() const { impl().print(); }

  Cons_Proxy cons() const { return impl().cons(); }
  Gens_Proxy gens() const { return impl().gens(); }
  Cons copy_cons() const { return impl().copy_cons(); }
  Gens copy_gens() const { return impl().copy_gens(); }
  Cons_Proxy normalized_cons() const { return impl().normalized_cons(); }
  Cons copy_normalized_cons() const { return impl().copy_normalized_cons(); }
  Gens_Info gens_info() const { return impl().gens_info(); }
  dim_type num_min_cons() const { return impl().num_min_cons(); }
  dim_type num_min_gens() const { return impl().num_min_gens(); }

  void collapse(dim_type n) { impl().collapse(n); }
  dim_type num_disjuncts() const { return impl().num_disjuncts(); }
  Cons_Proxy disjunct_cons(dim_type n) const {
    return impl().disjunct_cons(n);
  }
  bool geom_covers(const Poly& y) const {
    return impl().geom_covers(y.impl());
  }

  /* Modifiers */
  void m_swap(Poly& y) noexcept { impl().m_swap(y.impl()); }

  void set_empty() { impl().set_empty(); }
  void set_universe() { impl().set_universe(); }
  void set_topology(Topol t) { impl().set_topology(t); }
  void add_con(Con c) { impl().add_con(std::move(c)); }
  void add_cons(Cons cs) { impl().add_cons(std::move(cs)); }
  template <typename Iter>
  void add_cons(Iter first, Iter last) { impl().add_cons(first, last); }

  void add_gen(Gen g) { impl().add_gen(std::move(g)); }
  void add_gens(Gens gs) { impl().add_gens(std::move(gs)); }
  template <typename Iter>
  void add_gens(Iter first, Iter last) { impl().add_gens(first, last); }

  void topological_closure_assign() { impl().topological_closure_assign(); }

  template <typename Iter>
  void unconstrain(Iter first, Iter last) {
    impl().unconstrain(first, last);
  }
  void unconstrain(Var var) { impl().unconstrain(var); }
  void unconstrain(const Index_Set& vars) { impl().unconstrain(vars); }
  void intersection_assign(const Poly& y) {
    impl().intersection_assign(y.impl());
  }
  void join_assign(const Poly& y) { impl().join_assign(y.impl()); }
  void poly_hull_assign(const Poly& y) {
    impl().poly_hull_assign(y.impl());
  }
  void con_hull_assign(const Poly& y, bool boxed = false);
  void poly_difference_assign(const Poly& y);

  void affine_image(Var var, const Linear_Expr& expr,
                    const Integer& inhomo = Integer::zero(),
                    const Integer& den = Integer::one()) {
    impl().affine_image(var, expr, inhomo, den);
  }

  void affine_preimage(Var var, const Linear_Expr& expr,
                       const Integer& inhomo = Integer::zero(),
                       const Integer& den = Integer::one()) {
    impl().affine_preimage(var, expr, inhomo, den);
  }
  void
  parallel_affine_image(const Vars& vars,
                        const Linear_Exprs& exprs,
                        const Integers& inhomos, const Integers& dens) {
    impl().parallel_affine_image(vars, exprs, inhomos, dens);
  }

  void widening_assign(const Poly& y,
                       Widen_Impl w_impl = get_widen_impl(),
                       Widen_Spec w_spec = get_widen_spec()) {
    return impl().widening_assign(y.impl(), nullptr, w_impl, w_spec);
  }
  void widening_assign(const Poly& y, const Cons& upto_cons,
                       Widen_Impl w_impl = get_widen_impl(),
                       Widen_Spec w_spec = get_widen_spec()) {
    return impl().widening_assign(y.impl(), &upto_cons, w_impl, w_spec);
  }

  void time_elapse_assign(const Poly& y) {
    return impl().time_elapse_assign(y.impl());
  }

  /* Split */

  Poly split(const Con& c, Topol t) {
    Poly res;
    res.impl() = impl().split(c, t);
    return res;
  }
  Poly split(const Con& c) { return split(c, topology()); }

  Poly integral_split(const Con& c) {
    assert(is_necessarily_closed());
    Poly res;
    res.impl() = impl().integral_split(c);
    return res;
  }

  /* Change of space dim */

  void add_space_dims(dim_type m, bool project = false) {
    impl().add_space_dims(m, project);
  }
  void concatenate_assign(const Poly& y) {
    impl().concatenate_assign(y.impl());
  }
  void map_space_dims(const Dims& pfunc) {
    impl().map_space_dims(pfunc);
  }
  void remove_space_dim(Var var) { impl().remove_space_dim(var); }
  template <typename Iter>
  void remove_space_dims(Iter first, Iter last) {
    impl().remove_space_dims(first, last);
  }
  void remove_space_dims(const Index_Set& vars) {
    impl().remove_space_dims(vars);
  }
  void remove_higher_space_dims(dim_type new_dim) {
    impl().remove_higher_space_dims(new_dim);
  }
  void expand_space_dim(Var var, dim_type m) {
    impl().expand_space_dim(var, m);
  }
  void fold_space_dims(const Index_Set& vars, Var dest) {
    impl().fold_space_dims(vars, dest);
  }

  // This is semantically const,
  // but may affect syntactic representation.
  void minimize() const { impl().minimize(); }
  using Impl::get_minimize_filter_threshold;
  using Impl::set_minimize_filter_threshold;

  /* Input-output */
  bool ascii_load(std::istream& s) { return impl().ascii_load(s); }
  void ascii_dump(std::ostream& s) const { impl().ascii_dump(s); }

}; // class Poly

NOTHROW_MOVES(Poly_Impl);
NOTHROW_MOVES(Poly);

inline void swap(Poly& x, Poly& y) noexcept { x.m_swap(y); }

inline bool operator==(const Poly& x, const Poly& y) { return x.equals(y); }
inline bool operator!=(const Poly& x, const Poly& y) { return !(x == y); }

namespace IO_Operators {

std::ostream& operator<<(std::ostream& s, Topol topol);
std::ostream& operator<<(std::ostream& s, Poly_Impl::Status status);
std::ostream& operator<<(std::ostream& s, const Poly_Impl& ph);
std::ostream& operator<<(std::ostream& s, const Poly& ph);

} // namespace IO_Operators

template <typename Iter,
          typename PH = typename std::iterator_traits<Iter>::value_type>
void con_hull(PH& ph, Iter first, Iter last, bool boxed = false);

template <typename Iter,
          typename PH = typename std::iterator_traits<Iter>::value_type>
PH con_hull(Iter first, Iter last, bool boxed = false);

template <typename PH>
inline bool
boxed_is_disjoint_from(const PH& x, const PH& y,
                       const BBox& x_box, const BBox& y_box) {
  assert(x.topology() == y.topology());
  assert(x.space_dim() == y.space_dim());
  assert(x_box == x.get_bounding_box());
  assert(y_box == y.get_bounding_box());
  // Check boxes before.
  if (x_box.is_disjoint_from(y_box))
    return true;
  return x.is_disjoint_from(y);
}

} // namespace pplite

#endif // !defined(pplite_Poly_hh)
