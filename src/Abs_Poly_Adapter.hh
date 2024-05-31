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

#ifndef pplite_Abs_Poly_Adapter_hh
#define pplite_Abs_Poly_Adapter_hh 1

#include "Abs_Poly.hh"

#include "Poly.hh"
#include "B_Poly.hh"
#include "F_Poly.hh"
#include "U_Poly.hh"
#include "Poly_Stats.hh"
#include "PolySet.hh"
#include "memory_in_bytes.hh"

#include <cassert>
#include <iostream>
#include <string>
#include <utility>

namespace pplite {
namespace dynamic {

enum class Abs_Poly::Kind {
  POLY, POLY_STATS,
  B_POLY, B_POLY_STATS,
  F_POLY, F_POLY_STATS,
  U_POLY, U_POLY_STATS,
  UF_POLY, UF_POLY_STATS,
  P_SET, P_SET_STATS,
  FP_SET, FP_SET_STATS,
  SIZE
};

constexpr auto abs_poly_kind_size = static_cast<int>(Abs_Poly::Kind::SIZE);

// Template function mapping poly types to Kind.
template <typename PH>
constexpr Abs_Poly::Kind poly_type_to_kind();
// Explicit specializations.
template <> constexpr Abs_Poly::Kind
poly_type_to_kind<Poly>() { return Abs_Poly::Kind::POLY; }
template <> constexpr Abs_Poly::Kind
poly_type_to_kind<Stats<Poly>>() { return Abs_Poly::Kind::POLY_STATS; }
template <> constexpr Abs_Poly::Kind
poly_type_to_kind<B_Poly>() { return Abs_Poly::Kind::B_POLY; }
template <> constexpr Abs_Poly::Kind
poly_type_to_kind<Stats<B_Poly>>() { return Abs_Poly::Kind::B_POLY_STATS; }
template <> constexpr Abs_Poly::Kind
poly_type_to_kind<F_Poly>() { return Abs_Poly::Kind::F_POLY; }
template <> constexpr Abs_Poly::Kind
poly_type_to_kind<Stats<F_Poly>>() { return Abs_Poly::Kind::F_POLY_STATS; }
template <> constexpr Abs_Poly::Kind
poly_type_to_kind<U_Poly>() { return Abs_Poly::Kind::U_POLY; }
template <> constexpr Abs_Poly::Kind
poly_type_to_kind<Stats<U_Poly>>() { return Abs_Poly::Kind::U_POLY_STATS; }
template <> constexpr Abs_Poly::Kind
poly_type_to_kind<UF_Poly>() { return Abs_Poly::Kind::UF_POLY; }
template <> constexpr Abs_Poly::Kind
poly_type_to_kind<Stats<UF_Poly>>() { return Abs_Poly::Kind::UF_POLY_STATS; }
template <> constexpr Abs_Poly::Kind
poly_type_to_kind<P_Set>() { return Abs_Poly::Kind::P_SET; }
template <> constexpr Abs_Poly::Kind
poly_type_to_kind<Stats<P_Set>>() { return Abs_Poly::Kind::P_SET_STATS; }
template <> constexpr Abs_Poly::Kind
poly_type_to_kind<FP_Set>() { return Abs_Poly::Kind::FP_SET; }
template <> constexpr Abs_Poly::Kind
poly_type_to_kind<Stats<FP_Set>>() { return Abs_Poly::Kind::FP_SET_STATS; }


// Helpers mapping Kind to names (and viceversa).

constexpr const char* abs_poly_kind_names[abs_poly_kind_size] = {
  "Poly", "Poly_Stats",
  "B_Poly", "B_Poly_Stats",
  "F_Poly", "F_Poly_Stats",
  "U_Poly", "U_Poly_Stats",
  "UF_Poly", "UF_Poly_Stats",
  "P_Set", "P_Set_Stats",
  "FP_Set", "FP_Set_Stats"
};

constexpr const char*
abs_poly_kind_to_name(Abs_Poly::Kind kind) {
  // Let the header files be c++11 compliant.
  // Note: in c++11 constexpr functions cannot define local variables;
  // hence, we need to repeat the static cast three times.
  return (static_cast<int>(kind) >= 0 &&
          static_cast<int>(kind) < abs_poly_kind_size)
    ? abs_poly_kind_names[static_cast<int>(kind)]
    : "INVALID_KIND";
}

bool
abs_poly_name_to_kind(const char* cname,
                      Abs_Poly::Kind& kind,
                      bool& requires_stats);

// A template for a polymorphic sequence
template <typename Value, typename Sequence>
struct Abs_Sequence_Adapter : public Abs_Sequence<Value> {
private:
  using Base = Abs_Sequence<Value>;
  using seq_type = Sequence;
  seq_type seq;

public:
  using value_type = typename Base::value_type;

  explicit Abs_Sequence_Adapter(seq_type&& s) : seq(std::move(s)) {}
  ~Abs_Sequence_Adapter() override {}

  Abs_Sequence_Adapter(const Abs_Sequence_Adapter&) = delete;
  Abs_Sequence_Adapter(Abs_Sequence_Adapter&&) = delete;
  Abs_Sequence_Adapter& operator=(const Abs_Sequence_Adapter&) = delete;
  Abs_Sequence_Adapter& operator=(Abs_Sequence_Adapter&&) = delete;

  // Note: size() may be less than end_pos() due to skippable values.
  dim_type end_pos() const override { return seq.end_pos(); }
  dim_type size() const override { return seq.size(); }
  bool is_skippable(dim_type pos) const override {
    return seq.is_skippable(pos);
  }
  const value_type* get_value_ptr(dim_type pos) const override {
    return seq.get_value_ptr(pos);
  }
}; // Abs_Sequence_Adapter

template <typename PH>
class Abs_Poly_Adapter : public Abs_Poly {
public:
  using poly_type = PH;
  poly_type ph;

  // Helper
  static const poly_type& get_poly(const Abs_Poly& y) {
    assert(dynamic_cast<const Abs_Poly_Adapter*>(&y));
    return static_cast<const Abs_Poly_Adapter&>(y).ph;
  }

  explicit Abs_Poly_Adapter(const poly_type& y_ph) : ph(y_ph) {}
  explicit Abs_Poly_Adapter(poly_type&& y_ph) : ph(std::move(y_ph)) {}

  ~Abs_Poly_Adapter() = default;
  Abs_Poly_Adapter(const Abs_Poly_Adapter&) = delete;
  Abs_Poly_Adapter(Abs_Poly_Adapter&&) = delete;
  Abs_Poly_Adapter& operator=(const Abs_Poly_Adapter&) = delete;
  Abs_Poly_Adapter& operator=(Abs_Poly_Adapter&&) = delete;

  // Public interface
  explicit Abs_Poly_Adapter(dim_type d, Spec_Elem s, Topol t)
    : ph(d, s, t) {}

  Abs_Poly_Adapter* clone() const override {
    return new Abs_Poly_Adapter(ph);
  }

  void minimize() const override { ph.minimize(); }

  Kind poly_kind() const override {
    return poly_type_to_kind<poly_type>();
  }
  bool is_necessarily_closed() const override {
    return ph.is_necessarily_closed();
  }
  bool is_disjunctive() const override {
    return is_PolySet<poly_type>::value;
  }
  bool check_inv() const override { return ph.check_inv(); };

  bool is_empty() const override { return ph.is_empty(); }
  bool is_universe() const override { return ph.is_universe(); }
  bool is_bounded() const override { return ph.is_bounded(); }
  bool is_bounded_expr(bool from_below,
                       const Linear_Expr& expr) const override {
    return ph.is_bounded_expr(from_below, expr);
  }
  bool is_topologically_closed() const override {
    return ph.is_topologically_closed();
  }
  bool boxed_contains(const Abs_Poly& y) const override {
    return ph.boxed_contains(get_poly(y));
  }
  bool constrains(Var v) const override { return ph.constrains(v); }
  bool contains(const Abs_Poly& y) const override {
    return ph.contains(get_poly(y));
  }
  bool equals(const Abs_Poly& y) const override {
    return ph.equals(get_poly(y));
  }
  bool is_disjoint_from(const Abs_Poly& y) const override {
    return ph.is_disjoint_from(get_poly(y));
  }

  Topol topology() const override { return ph.topology(); }
  dim_type space_dim() const override { return ph.space_dim(); }
  dim_type affine_dim() const override { return ph.affine_dim(); }
  dim_type num_min_cons() const override { return ph.num_min_cons(); }
  dim_type num_min_gens() const override { return ph.num_min_gens(); }
  Cons copy_cons() const override { return ph.copy_cons(); }
  Gens copy_gens() const override { return ph.copy_gens(); }
  Poly_Con_Rel relation_with(const Con& c) const override {
    return ph.relation_with(c);
  }
  Poly_Gen_Rel relation_with(const Gen& g) const override {
    return ph.relation_with(g);
  }
  BBox get_bounding_box() const override { return ph.get_bounding_box(); }
  bool min(const Affine_Expr& ae, Rational& value,
           bool* included_ptr, Gen* g_ptr) const override {
    return ph.min(ae, value, included_ptr, g_ptr);
  }
  bool max(const Affine_Expr& ae, Rational& value,
           bool* included_ptr, Gen* g_ptr) const override {
    return ph.max(ae, value, included_ptr, g_ptr);
  }

  Itv get_bounds(Var var) const override { return ph.get_bounds(var); }
  Itv get_bounds(const Affine_Expr& ae) const override {
    return ph.get_bounds(ae);
  }
  Itv get_bounds(const Itv_Expr& ie) const override {
    return ph.get_bounds(ie);
  }

  Index_Set get_unconstrained() const override {
    return ph.get_unconstrained();
  }

  size_t hash() const override { return ph.hash(); }
  size_t get_memory_in_bytes() const override {
    return pplite::total_memory_in_bytes(ph);
  }

  // Sequence adapters
  using Con_Seq = Abs_Sequence_Adapter<Con, typename poly_type::Cons_Proxy>;
  using Gen_Seq = Abs_Sequence_Adapter<Gen, typename poly_type::Gens_Proxy>;

  // Proxies.
  Cons_Proxy cons() const override {
    return Cons_Proxy(new Con_Seq(ph.cons()));
  }
  Cons_Proxy normalized_cons() const override {
    return Cons_Proxy(new Con_Seq(ph.normalized_cons()));
  }
  Gens_Proxy gens() const override {
    return Gens_Proxy(new Gen_Seq(ph.gens()));
  }

  // Extension for powerset-like abstact domains.
  void collapse(dim_type n) override { return ph.collapse(n); }
  dim_type num_disjuncts() const override { return ph.num_disjuncts(); }
  Cons_Proxy disjunct_cons(dim_type n) const override {
    return Cons_Proxy(new Con_Seq(ph.disjunct_cons(n)));
  }
  bool geom_covers(const Abs_Poly& y) const override {
    return ph.geom_covers(get_poly(y));
  }

  // Modifiers
  void m_swap(Abs_Poly& y) noexcept override {
    assert(dynamic_cast<Abs_Poly_Adapter*>(&y));
    ph.m_swap(static_cast<Abs_Poly_Adapter&>(y).ph);
  }
  void set_empty() override { ph.set_empty(); }
  void set_universe() override { ph.set_universe(); }
  void set_topology(Topol t) override { ph.set_topology(t); }
  void add_con(const Con& c) override { ph.add_con(c); }
  void add_cons(const Cons& cs) override { ph.add_cons(cs); }
  void add_gen(const Gen& g) override { ph.add_gen(g); }
  void add_gens(const Gens& gs) override { ph.add_gens(gs); }

  void con_hull_assign(const Abs_Poly& y, bool boxed) override {
    ph.con_hull_assign(get_poly(y), boxed);
  }
  void concatenate_assign(const Abs_Poly& y) override {
    ph.concatenate_assign(get_poly(y));
  }
  void intersection_assign(const Abs_Poly& y) override {
    ph.intersection_assign(get_poly(y));
  }
  void join_assign(const Abs_Poly& y) override {
    ph.join_assign(get_poly(y));
  }
  void poly_hull_assign(const Abs_Poly& y) override {
    ph.poly_hull_assign(get_poly(y));
  }
  void poly_difference_assign(const Abs_Poly& y) override {
    ph.poly_difference_assign(get_poly(y));
  }
  void time_elapse_assign(const Abs_Poly& y) override {
    ph.time_elapse_assign(get_poly(y));
  }
  void topological_closure_assign() override {
    ph.topological_closure_assign();
  }
  void widening_assign(const Abs_Poly& y,
                       Widen_Impl w_impl, Widen_Spec w_spec) override {
    ph.widening_assign(get_poly(y), w_impl, w_spec);
  }
  void widening_assign(const Abs_Poly& y, const Cons& cs,
                       Widen_Impl w_impl, Widen_Spec w_spec) override {
    ph.widening_assign(get_poly(y), cs, w_impl, w_spec);
  }
  Abs_Poly_Adapter* split(const Con& con, Topol topol) override {
    return new Abs_Poly_Adapter(ph.split(con, topol));
  }
  Abs_Poly_Adapter* split(const Con& c) override {
    return new Abs_Poly_Adapter(ph.split(c, topology()));
  }
  Abs_Poly_Adapter* integral_split(const Con& con) override {
    return new Abs_Poly_Adapter(ph.integral_split(con));
  }

  void
  affine_image(Var var, const Linear_Expr& expr,
               const Integer& inhomo = Integer::zero(),
               const Integer& den = Integer::one()) override {
    ph.affine_image(var, expr, inhomo, den);
  }
  void
  affine_preimage(Var var, const Linear_Expr& expr,
                  const Integer& inhomo = Integer::zero(),
                  const Integer& den = Integer::one()) override {
    ph.affine_preimage(var, expr, inhomo, den);
  }
  void
  parallel_affine_image(const Vars& vars,
                        const Linear_Exprs& exprs,
                        const Integers& inhomos,
                        const Integers& dens) override {
    ph.parallel_affine_image(vars, exprs, inhomos, dens);
  }

  void unconstrain(Var var) override { ph.unconstrain(var); }
  void unconstrain(const Index_Set& vars) override {
    ph.unconstrain(vars.begin(), vars.end());
  }
  void add_space_dims(dim_type d, bool project) override {
    ph.add_space_dims(d, project);
  }
  void map_space_dims(const Dims& pfunc) override {
    ph.map_space_dims(pfunc);
  }
  void remove_space_dim(Var var) override { ph.remove_space_dim(var); }
  void remove_space_dims(const Index_Set& vars) override {
    ph.remove_space_dims(vars.begin(), vars.end());
  }
  void remove_higher_space_dims(dim_type new_dim) override {
    ph.remove_higher_space_dims(new_dim);
  }
  void expand_space_dim(Var var, dim_type m) override {
    ph.expand_space_dim(var, m);
  }
  void fold_space_dims(const Index_Set& vars, Var dest) override{
    ph.fold_space_dims(vars, dest);
  }

  /* Input-output */
  void print(std::ostream& os) const override { ph.print(os); }
  void ascii_dump(std::ostream& s) const override { ph.ascii_dump(s); }

}; // Abs_Poly_Adapter

// Declarations of explicit instantiations of the available concrete classes.
extern template class Abs_Poly_Adapter<Poly>;
extern template class Abs_Poly_Adapter<Stats<Poly>>;
extern template class Abs_Poly_Adapter<B_Poly>;
extern template class Abs_Poly_Adapter<Stats<B_Poly>>;
extern template class Abs_Poly_Adapter<F_Poly>;
extern template class Abs_Poly_Adapter<Stats<F_Poly>>;
extern template class Abs_Poly_Adapter<U_Poly>;
extern template class Abs_Poly_Adapter<Stats<U_Poly>>;
extern template class Abs_Poly_Adapter<UF_Poly>;
extern template class Abs_Poly_Adapter<Stats<UF_Poly>>;

extern template class Abs_Poly_Adapter<P_Set>;
extern template class Abs_Poly_Adapter<Stats<P_Set>>;
extern template class Abs_Poly_Adapter<FP_Set>;
extern template class Abs_Poly_Adapter<Stats<FP_Set>>;

// This is the "default" kind, used when not specifying it.
extern PPLITE_TLS Abs_Poly::Kind default_poly_kind;

inline const char* get_default_poly_kind_name() {
  return abs_poly_kind_to_name(default_poly_kind);
}

// Returns true (and sets global variable accordingly)
// if "name" is the name of a (supported) concrete class.
bool set_default_poly_kind(const char* name, bool noisy = false);

// A poor man's object factory.
// Use this to actually build concrete objects
inline Abs_Poly*
make_poly(Abs_Poly::Kind kind, dim_type d, Spec_Elem s, Topol t) {
  switch (kind) {
  case Abs_Poly::Kind::POLY:
    return new Abs_Poly_Adapter<Poly>(d, s, t);
  case Abs_Poly::Kind::POLY_STATS:
    return new Abs_Poly_Adapter<Stats<Poly>>(d, s, t);
  case Abs_Poly::Kind::B_POLY:
    return new Abs_Poly_Adapter<B_Poly>(d, s, t);
  case Abs_Poly::Kind::B_POLY_STATS:
    return new Abs_Poly_Adapter<Stats<B_Poly>>(d, s, t);
  case Abs_Poly::Kind::F_POLY:
    return new Abs_Poly_Adapter<F_Poly>(d, s, t);
  case Abs_Poly::Kind::F_POLY_STATS:
    return new Abs_Poly_Adapter<Stats<F_Poly>>(d, s, t);
  case Abs_Poly::Kind::U_POLY:
    return new Abs_Poly_Adapter<U_Poly>(d, s, t);
  case Abs_Poly::Kind::U_POLY_STATS:
    return new Abs_Poly_Adapter<Stats<U_Poly>>(d, s, t);
  case Abs_Poly::Kind::UF_POLY:
    return new Abs_Poly_Adapter<UF_Poly>(d, s, t);
  case Abs_Poly::Kind::UF_POLY_STATS:
    return new Abs_Poly_Adapter<Stats<UF_Poly>>(d, s, t);
  case Abs_Poly::Kind::P_SET:
    return new Abs_Poly_Adapter<P_Set>(d, s, t);
  case Abs_Poly::Kind::P_SET_STATS:
    return new Abs_Poly_Adapter<Stats<P_Set>>(d, s, t);
  case Abs_Poly::Kind::FP_SET:
    return new Abs_Poly_Adapter<FP_Set>(d, s, t);
  case Abs_Poly::Kind::FP_SET_STATS:
    return new Abs_Poly_Adapter<Stats<FP_Set>>(d, s, t);
  default:
    std::cerr << "Invalid Abs_Poly::Kind value" << std::endl;
    PPLITE_UNREACH;
  }
}

// Helper to use concrete kind names.
inline Abs_Poly*
make_poly(const char* kind_name, dim_type d, Spec_Elem s, Topol t) {
  auto kind = Abs_Poly::Kind::POLY;
  bool requires_stats = false;
  if (abs_poly_name_to_kind(kind_name, kind, requires_stats))
    return make_poly(kind, d, s, t);
  std::cerr << "Invalid Abs_Poly::Kind name" << std::endl;
  PPLITE_UNREACH;
  return nullptr;
}

} // namespace dynamic
} // namespace pplite

#endif // pplite_Abs_Poly_Adapter_hh
