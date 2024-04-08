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

#ifndef pplite_BBox_hh
#define pplite_BBox_hh 1

#include "globals.hh"
#include "Bits.hh"
#include "Gen.hh"
#include "Itv.hh"
#include "Rational.hh"

#include "Local_Stats.hh"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

namespace pplite {

// Topologically closed n-dim (rational) box.
template <bool KVI>
class Box {
public:
  static constexpr bool keep_volume_info = KVI;
  // For volume in Volume_Info:
  // volume.first = number of infinite bounds (i.e., rays)
  // volume.second = *pseudo* volume
  using Volume_Info = std::pair<dim_type, Rational>;

  Box(const Box& y) = default;
  Box(Box&& y) = default;
  Box& operator=(const Box& y) = default;
  Box& operator=(Box&& y) = default;
  ~Box() = default;

  // Either the empty or the universe box.
  explicit Box(dim_type sd, Spec_Elem se = Spec_Elem::UNIVERSE);

  // Construct from Box with/without volume info.
  explicit Box(const Box<not KVI>& y);

  // Volume info
  Volume_Info compute_volume_info() const;
  void maybe_update_volume_info() {
    if (keep_volume_info)
      volume = compute_volume_info();
  }
  dim_type num_rays() const { return volume.first; }
  const Rational& pseudo_volume() const {
    assert(is_bounded());
    return volume.second;
  }
  dim_type& num_rays() { return volume.first; }
  Rational& pseudo_volume() { return volume.second; }

  // Box methods: read only.
  bool check_inv() const;

  dim_type space_dim() const { return num_rows(itvs); }
  dim_type affine_dim() const;

  bool is_empty() const { return empty; }
  bool is_universe() const;
  bool is_bounded() const;
  bool is_bounded_expr(bool from_below, const Linear_Expr& expr) const;

  Topol topology() const { return Topol::CLOSED; }
  bool is_topologically_closed() const { return true; }

  bool constrains(Var v) const {
    return is_empty() || not itvs[v.id()].is_universe();
  }
  const Itv& get_bounds(Var v) const {
    assert(v.id() < space_dim());
    return is_empty() ? Itv::empty() : itvs[v.id()];
  }
  bool is_included_in(const Con& c) const;

  bool contains(const Box& y) const;
  bool equals(const Box& y) const;
  bool is_disjoint_from(const Box& y) const;

  // Note: x.less(y) implies !x.contains(y)
  bool less(const Box& y) const;

  dim_type num_min_cons() const;
  dim_type num_min_gens() const;
  Gens_Info gens_info() const;
  size_t hash() const;

  void print_itvs(std::ostream& os) const;
  void print(std::ostream& os) const;
  void ascii_dump(std::ostream& s) const;

  // Box methods: modifiers.
  void swap(Box& y) noexcept;
  void set_empty();
  void set_origin();
  void set_universe();
  void affine_image(Var var, Linear_Expr expr, Integer inhomo, Integer den);
  void concatenate_assign(const Box& y);
  void glb_assign(const Box& y);
  void lub_assign(const Box& y);
  void time_elapse_assign(const Box& y);
  void unconstrain(Var var);
  template <typename Iter>
  void unconstrain(Iter first, Iter last) {
    for ( ; first != last; ++first)
      unconstrain(Var(*first));
  }
  void unconstrain(const Index_Set& vars) {
    for (auto i : vars)
      unconstrain(Var(i));
  }
  void widening_assign(const Box& y);

  void add_space_dims(dim_type d, bool project = false);
  void permute_space_dims(const Dims& perm);
  void map_space_dims(const Dims& pfunc);
  void remove_space_dim(Var var) {
    remove_space_dims(Index_Set{var.id()});
  }
  void remove_space_dims(const Index_Set& vars);
  template <typename Iter>
  void remove_space_dims(Iter first, Iter last) {
    if (first == last)
      return;
    erase_using_sorted_indices(itvs, first, last);
    maybe_update_volume_info();
  }
  void remove_higher_space_dims(dim_type new_dim);
  void expand_space_dim(Var var, dim_type d);
  void fold_space_dims(const Index_Set& vars, Var dst);

  void add_gen(const Gen& g);
  template <typename Iter>
  void add_gens(Iter first, Iter last);
  void add_gens(const Gens& gs) { add_gens(gs.cbegin(), gs.cend()); }

  void refine_bounds(dim_type dim, const Itv& y_itv);
  void refine_as_integral(dim_type dim) {
    refine_as_integral(dim, dim+1);
  }
  void refine_as_integral(dim_type first, dim_type last);

  // Helpers: not really part of public interface
  void unconstrain_lb(dim_type d);
  void unconstrain_ub(dim_type d);

  // Modifier version of get_bounds().
  Itv& operator[](dim_type i) {
    assert(i < space_dim());
    return itvs[i];
  }

  // Read-only accessors for itvs.
  bool inf_lb(dim_type i) const { return itvs[i].inf_lb(); }
  bool inf_ub(dim_type i) const { return itvs[i].inf_ub(); }
  const Rational& lb(dim_type i) const { return itvs[i].lb; }
  const Rational& ub(dim_type i) const { return itvs[i].ub; }

public: // Meant, for the time being.
  bool empty;
  Itvs itvs;
  Volume_Info volume;
}; // class Box

template <bool KVI>
inline bool operator==(const Box<KVI>& x, const Box<KVI>& y) {
  return x.equals(y);
}

template <bool KVI>
inline bool operator!=(const Box<KVI>& x, const Box<KVI>& y) {
  return !(x == y);
}

template <bool KVI>
inline bool operator<(const Box<KVI>& x, const Box<KVI>& y) {
  return x.less(y);
}

template <bool KVI>
inline void
swap(Box<KVI>& x, Box<KVI>& y) noexcept {
  x.swap(y);
}

NOTHROW_MOVES(Box<true>);
NOTHROW_MOVES(Box<false>);

// BBox is a Box *without* volume info.
using BBox = Box<false>;

} // namespace pplite

#include "BBox_impl.hh"

#endif // !defined(pplite_BBox_hh)
