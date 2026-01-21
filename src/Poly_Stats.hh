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

#ifndef pplite_Poly_Stats_hh
#define pplite_Poly_Stats_hh 1

#include "globals.hh"
#include "clock.hh"
#include "Poly.hh"

namespace pplite {

// An enum providing a distinct value for each abstract operator.
// Note: listed in alphabetical order to avoid confusion;
// Z_OTHER_OPS stand for those operations that either cost nothing
// or are almost never called.
enum class AbsOp {
  ADD_CONS,
  ADD_GENS,
  AFFINE_IMAGES,
  AFFINE_IMAGES_PARALLEL,
  BOUNDS,
  BOXED_CONTAINS,
  CON_HULL,
  CONSTRAINS,
  CONTAINS,
  COPY_CONS,
  COPY_GENS,
  DIMS_ADD,
  DIMS_MAP,
  DIMS_REM,
  EQUALS,
  INTERSECTION,
  IS_DISJOINT,
  IS_EMPTY,
  JOIN,
  MINIMIZE,
  MINMAX,
  POLY_HULL,
  RELATION_WITH,
  SPLIT,
  TIME_ELAPSE,
  TOPCLOSURE,
  UNCONSTRAIN,
  WIDENING,
  WIDENING_UPTO,
  Z_OTHER_OPS,
  NUM_OPS
}; // enum class AbsOp

constexpr auto absop_size = static_cast<unsigned>(AbsOp::NUM_OPS);

// The array of abstract operator names (indexed by AbsOp).
// Careful: these must match the ordering of AbsOp enumeration values.
constexpr const char* absop_names[absop_size] = {
  "add_cons", "add_gens", "affine_images", "affine_images_parallel",
  "bounds", "boxed_contains",
  "con_hull", "constrains", "contains", "copy_cons", "copy_gens",
  "dims_add", "dims_map", "dims_rem", "equals",
  "intersection", "is_disjoint", "is_empty", "join",
  "minimize", "minmax", "poly_hull",
  "relation_with", "split", "time_elapse", "topclosure", "unconstrain",
  "widening", "widening_upto", "other_ops"
};

// A data structure collecting stats for each AbsOp.
struct AbsOp_Stats {
  // Whether or not the destructor should output statistics info.
  bool noisy_dtor;
  // The collected statistics.
  Time_Stats time_stats[absop_size];

  AbsOp_Stats() : noisy_dtor(false) {}

  // Clears all statistics.
  void reset();

  // Adds `time_incr' to stats for `op',
  // increasing by 1 the number of calls.
  void incr(AbsOp op, Clock::Duration time_incr) {
    time_stats[static_cast<unsigned>(op)].incr(time_incr);
  }

  void dump_op(std::ostream& os, AbsOp op, bool dump_times = false) const;
  void dump(std::ostream& os, bool dump_times = false) const;
  void dump_overall_time(std::ostream& os) const;
  ~AbsOp_Stats();
};

extern AbsOp_Stats absop_stats;

// If set to true, produces number of calls and computation time for
// polyhedra operators (when using the Stats<PH> wrapper!):
// the output of the stats occurs at end of program execution.
inline void set_noisy_stats(bool noisy) {
  absop_stats.noisy_dtor = noisy;
}

class AbsOp_Clock {
public:
  AbsOp_Clock(AbsOp op, bool noisy = false)
    : clock(), op_(op), noisy_(noisy) {}
  ~AbsOp_Clock() { stop_clock(); }
  // Neither copyable nor moveable.
  AbsOp_Clock(const AbsOp_Clock&) = delete;
  AbsOp_Clock& operator=(const AbsOp_Clock&) = delete;
  AbsOp_Clock(AbsOp_Clock&&) = delete;
  AbsOp_Clock& operator=(AbsOp_Clock&&) = delete;

  void stop_clock() const {
    auto elapsed = clock.elapsed_time();
    if (noisy_) {
      std::cerr << absop_names[static_cast<unsigned>(op_)] << " ";
      Clock::print(std::cerr, elapsed);
      std::cerr << std::endl;
    }
    absop_stats.incr(op_, elapsed);
  }

private:
  Clock clock;
  const AbsOp op_;
  const bool noisy_;
}; // class AbsOp_Clock

template <typename PH>
struct Stats : private PH {
  using Base = PH;

  using Base::check_inv;
  using Base::ascii_load;

  // Inherits the constructors: tracking accurate time stats
  // for these is quite tricky (and not really worth the effort).
  using Base::Base;

  Base& base() { return static_cast<Base&>(*this); }
  const Base& base() const { return static_cast<const Base&>(*this); }

  void add_con(Con c) {
    AbsOp_Clock clock(AbsOp::ADD_CONS);
    this->Base::add_con(std::move(c));
  }
  template <typename Iter>
  void add_cons(Iter first, Iter last) {
    AbsOp_Clock clock(AbsOp::ADD_CONS);
    this->Base::add_cons(first, last);
  }
  void add_cons(Cons cs) {
    AbsOp_Clock clock(AbsOp::ADD_CONS);
    this->Base::add_cons(std::move(cs));
  }
  void add_gen(Gen g) {
    AbsOp_Clock clock(AbsOp::ADD_GENS);
    this->Base::add_gen(std::move(g));
  }
  template <typename Iter>
  void add_gens(Iter first, Iter last) {
    AbsOp_Clock clock(AbsOp::ADD_GENS);
    this->Base::add_gens(first, last);
  }
  void add_gens(Gens gs) {
    AbsOp_Clock clock(AbsOp::ADD_GENS);
    this->Base::add_gens(std::move(gs));
  }
  void affine_image(Var var,
                    const Linear_Expr& expr,
                    const Integer& inhomo = Integer::zero(),
                    const Integer& den = Integer::one()) {
    AbsOp_Clock clock(AbsOp::AFFINE_IMAGES);
    this->Base::affine_image(var, expr, inhomo, den);
  }
  void affine_preimage(Var var,
                       const Linear_Expr& expr,
                       const Integer& inhomo = Integer::zero(),
                       const Integer& den = Integer::one()) {
    AbsOp_Clock clock(AbsOp::AFFINE_IMAGES);
    this->Base::affine_preimage(var, expr, inhomo, den);
  }
  void parallel_affine_image(const Vars& vars,
                             const Linear_Exprs& exprs,
                             const Integers& inhomos,
                             const Integers& dens) {
    AbsOp_Clock clock(AbsOp::AFFINE_IMAGES_PARALLEL);
    this->Base::parallel_affine_image(vars, exprs, inhomos, dens);
  }
  bool constrains(Var v) const {
    AbsOp_Clock clock(AbsOp::CONSTRAINS);
    return this->Base::constrains(v);
  }
  bool contains(const Stats& y) const {
    AbsOp_Clock clock(AbsOp::CONTAINS);
    return this->Base::contains(y);
  }
  bool strictly_contains(const Stats& y) const {
    AbsOp_Clock clock(AbsOp::CONTAINS);
    return this->Base::strictly_contains(y);
  }
  bool boxed_contains(const Stats& y) const {
    AbsOp_Clock clock(AbsOp::BOXED_CONTAINS);
    return this->Base::boxed_contains(y);
  }
  using Cons_Proxy = typename Base::Cons_Proxy;
  Cons_Proxy cons() const {
    return this->Base::cons();
  }
  using Gens_Proxy = typename Base::Gens_Proxy;
  Gens_Proxy gens() const {
    return this->Base::gens();
  }
  Cons copy_cons() const {
    AbsOp_Clock clock(AbsOp::COPY_CONS);
    return this->Base::copy_cons();
  }
  Gens copy_gens() const {
    AbsOp_Clock clock(AbsOp::COPY_GENS);
    return this->Base::copy_gens();
  }
  Cons_Proxy normalized_cons() const {
    return this->Base::normalized_cons();
  }
  Cons copy_normalized_cons() const {
    return this->Base::copy_normalized_cons();
  }
  template <typename Partial_Func>
  void map_space_dims(const Partial_Func& pfunc) {
    AbsOp_Clock clock(AbsOp::DIMS_MAP);
    this->Base::map_space_dims(pfunc);
  }
  void add_space_dims(dim_type m, bool project = false) {
    AbsOp_Clock clock(AbsOp::DIMS_ADD);
    this->Base::add_space_dims(m, project);
  }
  void remove_higher_space_dims(dim_type new_dim) {
    AbsOp_Clock clock(AbsOp::DIMS_REM);
    this->Base::remove_higher_space_dims(new_dim);
  }
  void remove_space_dim(Var v) {
    AbsOp_Clock clock(AbsOp::DIMS_REM);
    this->Base::remove_space_dim(v);
  }
  template <typename Iter>
  void remove_space_dims(Iter first, Iter last) {
    AbsOp_Clock clock(AbsOp::DIMS_REM);
    this->Base::remove_space_dims(first, last);
  }
  void remove_space_dims(const Index_Set& vs) {
    AbsOp_Clock clock(AbsOp::DIMS_REM);
    this->Base::remove_space_dims(vs);
  }
  bool equals(const Stats& y) const {
    AbsOp_Clock clock(AbsOp::EQUALS);
    return this->Base::equals(y);
  }
  BBox get_bounding_box() const {
    AbsOp_Clock clock(AbsOp::BOUNDS);
    return this->Base::get_bounding_box();
  }
  void intersection_assign(const Stats& y) {
    AbsOp_Clock clock(AbsOp::INTERSECTION);
    this->Base::intersection_assign(y);
  }
  bool is_disjoint_from(const Stats& y) const {
    AbsOp_Clock clock(AbsOp::IS_DISJOINT);
    return this->Base::is_disjoint_from(y);
  }
  bool is_empty() const {
    AbsOp_Clock clock(AbsOp::IS_EMPTY);
    return this->Base::is_empty();
  }
  void minimize() const {
    AbsOp_Clock clock(AbsOp::MINIMIZE);
    this->Base::minimize();
  }
  bool min(const Affine_Expr& ae, Rational& value,
           bool* included_ptr = nullptr, Gen* g_ptr = nullptr) const {
    AbsOp_Clock clock(AbsOp::MINMAX);
    return this->Base::min(ae, value, included_ptr, g_ptr);
  }
  bool max(const Affine_Expr& ae, Rational& value,
           bool* included_ptr = nullptr, Gen* g_ptr = nullptr) const {
    AbsOp_Clock clock(AbsOp::MINMAX);
    return this->Base::max(ae, value, included_ptr, g_ptr);
  }
  bool is_bounded_expr(bool from_below, const Linear_Expr& expr) const {
    AbsOp_Clock clock(AbsOp::BOUNDS);
    return this->Base::is_bounded_expr(from_below, expr);
  }
  Itv get_bounds(Var var) const {
    AbsOp_Clock clock(AbsOp::BOUNDS);
    return this->Base::get_bounds(var);
  }
  Itv get_bounds(const Affine_Expr& ae) const {
    AbsOp_Clock clock(AbsOp::BOUNDS);
    return this->Base::get_bounds(ae);
  }
  Itv get_bounds(const Itv_Expr& ie) const {
    AbsOp_Clock clock(AbsOp::BOUNDS);
    return this->Base::get_bounds(ie);
  }
  Index_Set get_unconstrained() const {
    AbsOp_Clock clock(AbsOp::CONSTRAINS);
    return this->Base::get_unconstrained();
  }
  void join_assign(const Stats& y) {
    AbsOp_Clock clock(AbsOp::JOIN);
    this->Base::join_assign(y);
  }
  void poly_hull_assign(const Stats& y) {
    AbsOp_Clock clock(AbsOp::POLY_HULL);
    this->Base::poly_hull_assign(y);
  }
  void con_hull_assign(const Stats& y, bool boxed = false) {
    AbsOp_Clock clock(AbsOp::CON_HULL);
    this->Base::con_hull_assign(y, boxed);
  }
  Poly_Con_Rel relation_with(const Con& c) const {
    AbsOp_Clock clock(AbsOp::RELATION_WITH);
    return this->Base::relation_with(c);
  }
  Poly_Gen_Rel relation_with(const Gen& g) const {
    AbsOp_Clock clock(AbsOp::RELATION_WITH);
    return this->Base::relation_with(g);
  }
  Stats split(const Con& c, Topol t) {
    AbsOp_Clock clock(AbsOp::SPLIT);
    Stats res;
    res.base() = this->Base::split(c, t);
    return res;
  }
  Stats split(const Con& c) { return split(c, topology()); }
  Stats integral_split(const Con& c) {
    AbsOp_Clock clock(AbsOp::SPLIT);
    Stats res;
    res.base() = this->Base::integral_split(c);
    return res;
  }
  void time_elapse_assign(const Stats& y) {
    AbsOp_Clock clock(AbsOp::TIME_ELAPSE);
    this->Base::time_elapse_assign(y);
  }
  void topological_closure_assign() {
    AbsOp_Clock clock(AbsOp::TOPCLOSURE);
    this->Base::topological_closure_assign();
  }
  void unconstrain(Var v) {
    AbsOp_Clock clock(AbsOp::UNCONSTRAIN);
    this->Base::unconstrain(v);
  }
  void unconstrain(const Index_Set& vs) {
    AbsOp_Clock clock(AbsOp::UNCONSTRAIN);
    this->Base::unconstrain(vs);
  }
  template <typename Iter>
  void unconstrain(Iter first, Iter last) {
    AbsOp_Clock clock(AbsOp::UNCONSTRAIN);
    this->Base::unconstrain(first, last);
  }
  void widening_assign(const Stats& y,
                       Widen_Impl w_impl = get_widen_impl(),
                       Widen_Spec w_spec = get_widen_spec()) {
    AbsOp_Clock clock(AbsOp::WIDENING);
    this->Base::widening_assign(y, w_impl, w_spec);
  }
  void widening_assign(const Stats& y, const Cons& upto_cons,
                       Widen_Impl w_impl = get_widen_impl(),
                       Widen_Spec w_spec = get_widen_spec()) {
    AbsOp_Clock clock(AbsOp::WIDENING_UPTO);
    this->Base::widening_assign(y, upto_cons, w_impl, w_spec);
  }
  // Other ops.
  dim_type affine_dim() const {
    AbsOp_Clock clock(AbsOp::Z_OTHER_OPS);
    return this->Base::affine_dim();
  }
  void concatenate_assign(const Stats& y) {
    AbsOp_Clock clock(AbsOp::Z_OTHER_OPS);
    this->Base::concatenate_assign(y);
  }
  void expand_space_dim(Var v, dim_type n) {
    AbsOp_Clock clock(AbsOp::Z_OTHER_OPS);
    this->Base::expand_space_dim(v, n);
  }
  void fold_space_dims(const Index_Set& vs, Var v) {
    AbsOp_Clock clock(AbsOp::Z_OTHER_OPS);
    this->Base::fold_space_dims(vs, v);
  }
  size_t hash() const {
    AbsOp_Clock clock(AbsOp::Z_OTHER_OPS);
    return this->Base::hash();
  }
  void m_swap(Stats& y) noexcept {
    AbsOp_Clock clock(AbsOp::Z_OTHER_OPS);
    this->Base::m_swap(y);
  }
  void set_empty() {
    AbsOp_Clock clock(AbsOp::Z_OTHER_OPS);
    this->Base::set_empty();
  }
  void set_universe() {
    AbsOp_Clock clock(AbsOp::Z_OTHER_OPS);
    this->Base::set_universe();
  }
  bool is_bounded() const {
    AbsOp_Clock clock(AbsOp::Z_OTHER_OPS);
    return this->Base::is_bounded();
  }
  bool is_topologically_closed() const {
    AbsOp_Clock clock(AbsOp::Z_OTHER_OPS);
    return this->Base::is_topologically_closed();
  }
  bool is_necessarily_closed() const {
    AbsOp_Clock clock(AbsOp::Z_OTHER_OPS);
    return this->Base::is_necessarily_closed();
  }
  bool is_universe() const {
    AbsOp_Clock clock(AbsOp::Z_OTHER_OPS);
    return this->Base::is_universe();
  }
  dim_type num_min_cons() const {
    AbsOp_Clock clock(AbsOp::Z_OTHER_OPS);
    return this->Base::num_min_cons();
  }
  dim_type num_min_gens() const {
    AbsOp_Clock clock(AbsOp::Z_OTHER_OPS);
    return this->Base::num_min_gens();
  }
  void poly_difference_assign(const Stats& y) {
    AbsOp_Clock clock(AbsOp::Z_OTHER_OPS);
    this->Base::poly_difference_assign(y);
  }
  dim_type space_dim() const {
    AbsOp_Clock clock(AbsOp::Z_OTHER_OPS);
    return this->Base::space_dim();
  }
  Topol topology() const {
    AbsOp_Clock clock(AbsOp::Z_OTHER_OPS);
    return this->Base::topology();
  }
  void set_topology(Topol t) {
    AbsOp_Clock clock(AbsOp::Z_OTHER_OPS);
    this->Base::set_topology(t);
  }
  size_t total_memory_in_bytes() const {
    AbsOp_Clock clock(AbsOp::Z_OTHER_OPS);
    return this->Base::total_memory_in_bytes();
  }
  void print(std::ostream& os) const {
    AbsOp_Clock clock(AbsOp::Z_OTHER_OPS);
    this->Base::print(os);
  }
  void print() const {
    AbsOp_Clock clock(AbsOp::Z_OTHER_OPS);
    this->Base::print();
  }
  void ascii_dump(std::ostream& os) const {
    this->Base::ascii_dump(os);
  }

  void collapse(dim_type n) {
    this->Base::collapse(n);
  }
  dim_type num_disjuncts() const {
    return this->Base::num_disjuncts();
  }
  Cons_Proxy disjunct_cons(dim_type n) const {
    return this->Base::disjunct_cons(n);
  }
  bool geom_covers(const Stats& y) const {
    AbsOp_Clock clock(AbsOp::CONTAINS);
    return this->Base::geom_covers(y);
  }

}; // struct Stats

template <typename PH>
bool operator==(const Stats<PH>& x, const Stats<PH>& y) {
  AbsOp_Clock clock(AbsOp::EQUALS);
  const PH& xx = x.base();
  const PH& yy = y.base();
  return xx == yy;
}

template <typename PH>
void swap(Stats<PH>& x, Stats<PH>& y) noexcept { x.m_swap(y); }

template <typename PH>
inline bool
operator!=(const Stats<PH>& x, const Stats<PH>& y) {
  return !(x == y);
}

namespace IO_Operators {

template <typename PH>
std::ostream& operator<<(std::ostream& s, const Stats<PH>& ph) {
  return s << ph.base();
}

} // namespace IO_Operators

template <typename PH>
inline bool
boxed_is_disjoint_from(const Stats<PH>& x, const Stats<PH>& y,
                       const BBox& x_box, const BBox& y_box) {
  AbsOp_Clock clock(AbsOp::IS_DISJOINT);
  const PH& xx = x.base();
  const PH& yy = y.base();
  return boxed_is_disjoint_from(xx, yy, x_box, y_box);
}

using Poly_Stats = Stats<Poly>;
NOTHROW_MOVES(Poly_Stats);

} // namespace pplite

#endif // !defined(pplite_Poly_Stats_hh)
