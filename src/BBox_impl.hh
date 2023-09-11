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


#ifndef pplite_BBox_impl_hh
#define pplite_BBox_impl_hh 1

namespace pplite {

template <bool KVI>
inline
Box<KVI>::Box(dim_type sd, Spec_Elem se)
  : empty(se == Spec_Elem::EMPTY),
    itvs(sd),
    volume(0, Rational::zero()) {
  if (keep_volume_info) {
    num_rays() = empty ? 0 : 2 * space_dim();
    pseudo_volume() = (empty || space_dim() != 0)
      ? Rational::zero() : Rational::one();
  }
  assert(check_inv());
}

template <bool KVI>
inline
Box<KVI>::Box(const Box<not KVI>& y)
  : empty(y.empty),
    itvs(y.itvs),
    volume(0, Rational::zero()) {
  maybe_update_volume_info();
  assert(check_inv());
}

template <bool KVI>
inline bool
Box<KVI>::is_universe() const {
  if (keep_volume_info)
    return num_rays() == 2 * space_dim();
  return not is_empty() && all_of(itvs, std::mem_fn(&Itv::is_universe));
}

template <bool KVI>
inline bool
Box<KVI>::is_bounded() const {
  if (keep_volume_info)
    return num_rays() == 0;
  return is_empty() || all_of(itvs, std::mem_fn(&Itv::is_bounded));
}

template <bool KVI>
bool
Box<KVI>::is_bounded_expr(bool from_below, const Linear_Expr& expr) const {
  if (keep_volume_info && num_rays() == 0)
    return true;
  if (is_empty())
    return true;
  for (auto i : bwd_dim_range(expr)) {
    auto s = sgn(expr.get(i));
    if (s == 0)
      continue;
    if (s > 0) {
      if (inf_lb(i) && !from_below)
        return false;
      if (inf_ub(i) && from_below)
        return false;
    } else {
      assert(s < 0);
      if (inf_lb(i) && from_below)
        return false;
      if (inf_ub(i) && !from_below)
        return false;
    }
  }
  return true;
}

template <bool KVI>
bool
Box<KVI>::is_included_in(const Con& c) const {
  // If `c' is satisfied by the closest point, the whole box is.
  const auto& c_ex = c.linear_expr();
  Rational scal_prod(c.inhomo_term());
  for (auto i : bwd_dim_range(c)) {
    const auto& coeff = c_ex.get(i);
    const auto s = sgn(coeff);
    if (s == 0)
      continue;
    const auto& itv = itvs[i];
    // Check `itv' infinite bounds.
    if ((s > 0 && itv.inf_lb()) || (s < 0 && itv.inf_ub())) {
      // Coeff should have been 0.
      return false;
    }
    if (c.is_equality() && !itv.is_singleton())
      return false;
    // Keep building the scalar product.
    add_mul_assign(scal_prod, ((s > 0) ? itv.lb : itv.ub), Rational(coeff));
  }
  const auto sp_sign = sgn(scal_prod);
  if (c.is_equality())
    return sp_sign == 0;
  else if (c.is_strict_inequality())
    return sp_sign > 0;
  else
    return sp_sign >= 0;
}

template <bool KVI>
inline bool
Box<KVI>::contains(const Box& y) const {
  auto& x = *this;
  assert(x.space_dim() == y.space_dim());
  if (y.is_empty())
    return true;
  if (x.is_empty())
    return false;
  if (keep_volume_info) {
    if (x.is_bounded() && y.is_bounded()) {
      if (x.pseudo_volume() < y.pseudo_volume())
        return false;
    } else if (x.num_rays() < y.num_rays())
      return false;
  }
  // Note: uses std::equal, but actually checks for containment.
  return std::equal(x.itvs.begin(), x.itvs.end(),
                    y.itvs.begin(),
                    std::mem_fn<bool(const Itv&) const>(&Itv::contains));
}

template <bool KVI>
inline bool
Box<KVI>::equals(const Box& y) const {
  const auto& x = *this;
  if (x.space_dim() != y.space_dim())
    return false;
  if (x.is_empty() != y.is_empty())
    return false;
  if (x.is_empty())
    return true;
  if (keep_volume_info) {
    if (x.is_bounded() && y.is_bounded()) {
      if (x.pseudo_volume() != y.pseudo_volume())
        return false;
    } else if (x.num_rays() != y.num_rays())
      return false;
  }
  return (x.itvs == y.itvs);
}

template <bool KVI>
inline bool
Box<KVI>::less(const Box& y) const {
  const auto& x = *this;
  if (keep_volume_info)
    return x.volume < y.volume;
  if (x.is_bounded()) {
    if (!y.is_bounded())
      return true;
    switch (compare(x.pseudo_volume(), y.pseudo_volume())) {
    case -1:
      return true;
    case 1:
      return false;
    case 0:
    default:
      return x.itvs < y.itvs;
    }
  }
  assert(!x.is_bounded());
  if (y.is_bounded())
    return false;
  return x.itvs < y.itvs;
}

template <bool KVI>
inline bool
Box<KVI>::is_disjoint_from(const Box& y) const {
  const auto& x = *this;
  assert(x.space_dim() == y.space_dim());
  return x.is_empty() || y.is_empty()
    || not std::equal(x.itvs.begin(), x.itvs.end(),
                      y.itvs.begin(), std::mem_fn(&Itv::intersects));
}

template <bool KVI>
inline dim_type
Box<KVI>::affine_dim() const {
  if (is_empty())
    return 0;
  auto num_eqns = std::count_if(itvs.begin(), itvs.end(),
                                std::mem_fn(&Itv::is_singleton));
  return space_dim() - num_eqns;
}

template <bool KVI>
inline dim_type
Box<KVI>::num_min_cons() const {
  if (is_empty())
    return 1;
  auto res = 0;
  for (const auto& itv : itvs)
    res += itv.num_min_cons();
  return res;
}

template <bool KVI>
inline Gens_Info
Box<KVI>::gens_info() const {
  if (is_empty())
    return { 0, 0, 0, 0, 0 };
  auto ln = 0, r = 0, skp = 1;
  for (const auto& itv : itvs) {
    assert(not itv.is_empty());
    if (itv.is_universe())
      ++ln;
    else if (not itv.is_bounded())
      ++r;
    else if (not itv.is_singleton())
      skp *= 2;
  }
  return { ln, r, 0, skp, 0 };
}

template <bool KVI>
inline dim_type
Box<KVI>::num_min_gens() const {
  if (is_empty())
    return 0;
  auto gi = gens_info();
  return std::get<0>(gi) + std::get<1>(gi) + std::get<3>(gi);
}

template <bool KVI>
inline size_t
Box<KVI>::hash() const {
  auto res = hash_size(space_dim());
  if (is_empty())
    return res;
  for (const auto& itv : itvs)
    hash_combine(res, itv.hash());
  return res;
}

template <bool KVI>
inline void
Box<KVI>::print_itvs(std::ostream& os) const {
  for (auto i : index_range(itvs)) {
    if (i > 0)
      os << ", ";
    Var(i).print(os);
    os << " in ";
    itvs[i].print(os);
  }
}

template <bool KVI>
inline void
Box<KVI>::print(std::ostream& os) const {
  if (is_empty())
    os << "empty";
  else
    print_itvs(os);
}

template <bool KVI>
inline void
Box<KVI>::swap(Box& y) noexcept {
  auto& x =*this;
  using std::swap;
  swap(x.empty, y.empty);
  swap(x.itvs, y.itvs);
  if (keep_volume_info)
    swap(x.volume, y.volume);
}

template <bool KVI>
inline void
Box<KVI>::set_empty() {
  empty = true;
  if (keep_volume_info) {
    num_rays() = 0;
    pseudo_volume() = Rational::zero();
  }
  assert(check_inv());
}

template <bool KVI>
inline void
Box<KVI>::set_origin() {
  empty = false;
  for (auto& itv : itvs)
    itv.set_zero();
  if (keep_volume_info) {
    num_rays() = 0;
    pseudo_volume() = Rational::one();
  }
  assert(check_inv());
}

template <bool KVI>
inline void
Box<KVI>::set_universe() {
  empty = false;
  for (auto& itv : itvs)
    itv.set_universe();
  if (keep_volume_info) {
    num_rays() = 2 * space_dim();
    pseudo_volume() = (space_dim() == 0)
      ? Rational::one() : Rational::zero();
  }
}

template <bool KVI>
inline void
Box<KVI>::concatenate_assign(const Box& y) {
  auto& x = *this;
  assert(x.check_inv());
  assert(y.check_inv());
  x.itvs.reserve(x.space_dim() + y.space_dim());
  x.itvs.insert(x.itvs.end(), y.itvs.begin(), y.itvs.end());
  if (y.is_empty())
    x.set_empty();
  else if (keep_volume_info) {
    x.num_rays() += y.num_rays();
    if (x.num_rays() == 0)
      x.pseudo_volume() *= y.pseudo_volume();
    else
      x.pseudo_volume() = Rational::zero();
  }
  assert(x.check_inv());
}

template <bool KVI>
inline void
Box<KVI>::glb_assign(const Box& y) {
  auto& x = *this;
  assert(x.space_dim() == y.space_dim());
  assert(x.check_inv());
  assert(y.check_inv());
  if (x.is_empty())
    return;
  if (y.is_empty()) {
    x.set_empty();
    return;
  }
  for (auto i : bwd_dim_range(x)) {
    if (x.itvs[i].glb_assign(y.itvs[i])) {
      x.set_empty();
      return;
    }
  }
  x.maybe_update_volume_info();
  assert(x.check_inv());
}

template <bool KVI>
inline void
Box<KVI>::lub_assign(const Box& y) {
  auto& x = *this;
  assert(x.space_dim() == y.space_dim());
  assert(x.check_inv());
  assert(y.check_inv());
  if (y.is_empty())
    return;
  if (x.is_empty()) {
    x = y;
    return;
  }
  for (auto i : bwd_dim_range(x))
    x.itvs[i].lub_assign(y.itvs[i]);
  x.maybe_update_volume_info();
  assert(!x.is_empty());
  assert(x.check_inv());
}

template <bool KVI>
inline void
Box<KVI>::time_elapse_assign(const Box& y) {
  auto& x = *this;
  assert(x.space_dim() == y.space_dim());
  assert(x.check_inv());
  assert(y.check_inv());
  if (x.is_empty())
    return;
  if (y.is_empty()) {
    x.set_empty();
    return;
  }
  for (auto i : bwd_dim_range(x)) {
    auto& x_i = x.itvs[i];
    const auto& y_i = y.itvs[i];
    if (x_i.has_lb() && (y_i.inf_lb() || sgn(y_i.lb) < 0)) {
      x_i.unset_lb();
      if (keep_volume_info)
        ++x.num_rays();
    }
    if (x_i.has_ub() && (y_i.inf_ub() || sgn(y_i.ub) > 0)) {
      x_i.unset_ub();
      if (keep_volume_info)
        ++x.num_rays();
    }
  }
  if (keep_volume_info) {
    if (num_rays() > 0)
      pseudo_volume() = Rational::zero();
  }
  assert(!x.is_empty());
  assert(x.check_inv());
}

template <bool KVI>
inline void
Box<KVI>::unconstrain(Var var) {
  if (is_empty())
    return;
  auto& itv = itvs[var.id()];
  if (itv.is_universe())
    return;
  if (keep_volume_info) {
    num_rays() += (itv.is_bounded() ? 2 : 1);
    pseudo_volume() = Rational::zero();
  }
  itv.set_universe();
}

template <bool KVI>
inline void
Box<KVI>::widening_assign(const Box& y) {
  auto& x = *this;
  assert(x.check_inv() && y.check_inv());
  assert(x.space_dim() == y.space_dim());
  // Adopt safe implementation (no inclusion precondition).
  if (y.is_empty())
    return;
  if (x.is_empty()) {
    x = y;
    return;
  }
  if (x.space_dim() == 0) {
    assert(x.is_universe() && y.is_universe());
    return;
  }
  for (auto i : bwd_dim_range(x))
    x.itvs[i].widen_assign(y.itvs[i]);
  x.maybe_update_volume_info();
  assert(x.check_inv());
}

template <bool KVI>
inline void
Box<KVI>::add_space_dims(dim_type d, bool project) {
  if (d == 0)
    return;
  auto old_sd = space_dim();
  auto new_sd = old_sd + d;
  itvs.resize(new_sd);
  if (is_empty())
    return;
  if (project) {
    for (auto i : range(old_sd, new_sd))
      itvs[i].set_zero();
  } else {
    if (keep_volume_info) {
      num_rays() += 2*d;
      pseudo_volume() = Rational::zero();
    }
  }
  assert(check_inv());
}

template <bool KVI>
inline void
Box<KVI>::permute_space_dims(const Dims& perm) {
  assert(num_rows(perm) == space_dim());
  auto sd = space_dim();
  if (is_empty() || sd < 2)
    return;
  Itvs perm_itvs(sd);
  for (auto i : range(sd))
    perm_itvs[perm[i]] = std::move(itvs[i]);
  using std::swap;
  swap(itvs, perm_itvs);
  assert(check_inv());
}

template <bool KVI>
void
Box<KVI>::map_space_dims(const Dims& pfunc) {
  assert(num_rows(pfunc) == space_dim());
  const auto old_dim = space_dim();
  if (old_dim == 0)
    return;
  const auto max_dim = *std::max_element(pfunc.begin(), pfunc.end());
  const auto new_dim = (max_dim == not_a_dim()) ? 0 : (1 + max_dim);
  if (new_dim == 0 || is_empty()) {
    remove_higher_space_dims(new_dim);
    assert(check_inv());
    return;
  }

  if (new_dim == old_dim) {
    // pfunc is a permutation
    permute_space_dims(pfunc);
    assert(check_inv());
    return;
  }

  // `pfunc' is not a permutation: remove_space_dims + permutation.
  Index_Set tbr;
  Dims perm;
  for (auto i : range(old_dim)) {
    if (pfunc[i] == not_a_dim())
      tbr.set(i);
    else
      perm.push_back(pfunc[i]);
  }
  remove_space_dims(tbr);
  permute_space_dims(perm);
  assert(check_inv());
}

template <bool KVI>
inline void
Box<KVI>::remove_space_dims(const Index_Set& vars) {
  if (vars.empty())
    return;
  assert(vars.last() < space_dim());
  auto new_dim = space_dim() - vars.size();
  if (is_empty()) {
    itvs.resize(new_dim);
    assert(check_inv());
    return;
  }
  erase_using_sorted_indices(itvs, vars);
  maybe_update_volume_info();
  assert(check_inv());
}

template <bool KVI>
inline void
Box<KVI>::remove_higher_space_dims(dim_type new_dim) {
  assert(new_dim <= space_dim());
  if (new_dim == space_dim())
    return;
  itvs.resize(new_dim);
  if (!is_empty())
    maybe_update_volume_info();
  assert(check_inv());
}

template <bool KVI>
inline void
Box<KVI>::expand_space_dim(Var var, dim_type d) {
  assert(var.id() < space_dim());
  if (d == 0)
    return;
  if (is_empty()) {
    add_space_dims(d);
    return;
  }
  const auto& itv = itvs[var.id()];
  if (keep_volume_info) {
    if (num_rays() == 0) {
      auto itv_len = itv.length();
      if (!itv_len.is_zero()) {
        itv_len += Rational::one();
        pseudo_volume() *= pow_si(itv_len, d);
      }
    } else {
      num_rays() += d * itv.num_rays();
    }
  }
  itvs.resize(space_dim() + d, itv);
  assert(check_inv());
}

template <bool KVI>
inline void
Box<KVI>::fold_space_dims(const Index_Set& vars, Var dst) {
  assert(dst.id() < space_dim());
  assert(!vars.test(dst.id()));
  if (vars.empty())
    return;
  assert(vars.last() < space_dim());
  if (is_empty()) {
    itvs.resize(space_dim() - vars.size());
    return;
  }
  auto& dst_itv = itvs[dst.id()];
  for (auto i : vars)
    dst_itv.lub_assign(itvs[i]);
  erase_using_sorted_indices(itvs, vars);
  maybe_update_volume_info();
  assert(check_inv());
}

template <bool KVI>
inline void
Box<KVI>::ascii_dump(std::ostream& os) const {
  os << "empty = " << empty << ", ";
  os << "space_dim = " << space_dim() << "\n";
  if (keep_volume_info) {
    os << "num_inf = " << volume.first << ", ";
    os << "pseudo_vol = " << volume.second << "\n";
  }
  print_itvs(os);
  os << "\n";
}

template <bool KVI>
typename Box<KVI>::Volume_Info
Box<KVI>::compute_volume_info() const {
  Volume_Info res = { 0, Rational::zero() };
  if (is_empty())
    return res;
  for (const auto& itv : itvs)
    res.first += itv.num_rays();
  if (res.first == 0) {
    res.second = Rational::one();
    for (const auto& itv : itvs)
      res.second *= (Rational::one() + itv.length());
  }
  return res;
}

template <bool KVI>
bool
Box<KVI>::check_inv() const {
  if (keep_volume_info) {
    if (space_dim() == 0) {
      // Note: 0-dim bbox is bounded even when is UNIVERSE.
      if (!is_bounded())
        return false;
      if (empty)
        return pseudo_volume() == Rational::zero();
      else
        return pseudo_volume() == Rational::one();
    }
    if (!is_bounded() && volume.second != Rational::zero())
      return false;
  }

  if (is_empty()) {
    // itvs can be anything for empty boxes.
    return true;
  }
  if (not all_of(itvs, std::mem_fn(&Itv::check_inv)))
    return false;
  if (keep_volume_info) {
    if (volume != compute_volume_info())
      return false;
  }
  return true;
}

template <bool KVI>
inline void
Box<KVI>::unconstrain_lb(dim_type d) {
  assert(d < space_dim());
  auto& itv = itvs[d];
  if (!itv.has_lb())
    return;
  itv.unset_lb();
  if (keep_volume_info) {
    ++num_rays();
    pseudo_volume() = Rational::zero();
  }
}

template <bool KVI>
inline void
Box<KVI>::unconstrain_ub(dim_type d) {
  assert(d < space_dim());
  auto& itv = itvs[d];
  if (!itv.has_ub())
    return;
  itv.unset_ub();
  if (keep_volume_info) {
    ++num_rays();
    pseudo_volume() = Rational::zero();
  }
}

namespace detail {

template <bool KVI>
inline void
add_as_line(Box<KVI>& box, const Gen& g) {
  assert(g.space_dim() <= box.space_dim());
  assert(!box.is_empty());
  for (auto i : bwd_dim_range(g)) {
    if (g.coeff(Var(i)) == 0)
      continue;
    box.unconstrain(Var(i));
  }
}

template <bool KVI>
inline void
add_as_ray(Box<KVI>& box, const Gen& g) {
  assert(g.space_dim() <= box.space_dim());
  assert(!box.is_empty());
  for (auto i : bwd_dim_range(g)) {
    auto s = sgn(g.coeff(Var(i)));
    if (s < 0)
      box.unconstrain_lb(i);
    else if (s > 0)
      box.unconstrain_ub(i);
  }
}

template <bool KVI>
inline void
add_as_point(Box<KVI>& box, const Gen& g) {
  assert(g.space_dim() <= box.space_dim());
  assert(!box.is_empty());

  for (auto i : bwd_dim_range(box)) {
    if (not box.constrains(Var(i)))
      continue;
    Rational value(g.coeff(Var(i)), g.divisor());
    auto& itv = box[i];
    if (itv.has_lb() && value < itv.lb)
      itv.lb = value;
    if (itv.has_ub() && value > itv.ub)
      itv.ub = value;
  }
}

template <bool KVI>
void
add_gen_unsafe(Box<KVI>& box, const Gen& g) {
  // Note: assumes box not empty and does not update volume info,
  // so that invariant may hot hold.
  assert(!box.is_empty());
  assert(g.space_dim() <= box.space_dim());
  assert(!(box.space_dim() == 0 && g.is_line_or_ray()));

  if (box.space_dim() == 0)
    return;

  auto type = g.type();
  switch (type) {
  case Gen::LINE:
    detail::add_as_line(box, g);
    break;
  case Gen::RAY:
    detail::add_as_ray(box, g);
    break;
  case Gen::POINT:
  case Gen::CLOSURE_POINT:
  default:
    assert(type == Gen::POINT || type == Gen::CLOSURE_POINT);
    detail::add_as_point(box, g);
    break;
  }
}

template <bool KVI>
inline void
init_with_point(Box<KVI>& box, const Gen& g) {
  assert(g.space_dim() <= box.space_dim());
  assert(g.is_point());
  assert(box.is_empty());

  box.set_origin();
  for (auto i : bwd_dim_range(g)) {
    const auto& c_i = g.coeff(Var(i));
    if (c_i.is_zero())
      continue;
    box[i].set_singleton(Rational(c_i, g.divisor()));
  }
  if (box.keep_volume_info) {
    box.num_rays() = 0;
    box.pseudo_volume() = Rational::one();
  }
}

} // namespace detail

template <bool KVI>
void
Box<KVI>::add_gen(const Gen& g) {
  assert(g.space_dim() <= space_dim());
  assert(!(space_dim() == 0 && g.is_line_or_ray()));

  if (space_dim() == 0) {
    if (is_empty()) {
      assert(g.is_point());
      empty = false;
      pseudo_volume() = Rational::one();
    }
    assert(check_inv());
    return;
  }

  auto& x = *this;
  if (is_empty())
    detail::init_with_point(x, g);
  else {
    detail::add_gen_unsafe(x, g);
    maybe_update_volume_info();
  }
  assert(!is_empty());
  assert(check_inv());
}

template <bool KVI>
template <typename Iter>
void
Box<KVI>::add_gens(Iter first, Iter last) {
  auto& x = *this;
  if (first == last)
    return;
  if (x.is_empty()) {
    auto iter = std::find_if(first, last,
                             std::mem_fn(&Gen::is_point));
    assert(iter != last);
    x.add_gen(*iter);
    x.add_gens(first, iter);
    ++iter;
    x.add_gens(iter, last);
    return;
  }

  assert(!x.is_empty());
  for ( ; first != last; ++first)
    detail::add_gen_unsafe(x, *first);
  x.maybe_update_volume_info();
  assert(check_inv());
}

template <bool KVI>
void
Box<KVI>::refine_bounds(dim_type dim, const Itv& y_itv) {
  assert(check_inv() && !is_empty());
  assert(dim < space_dim());
  auto& x_itv = itvs[dim];
  if (y_itv.contains(x_itv))
    return;
  if (x_itv.glb_assign(y_itv))
    set_empty();
  else
    maybe_update_volume_info();
  assert(check_inv());
}

template <bool KVI>
void
Box<KVI>::refine_as_integral(dim_type first, dim_type last) {
  auto& x = *this;
  assert(x.check_inv());
  assert(0 <= first && first <= last && last <= x.space_dim());
  if (x.is_empty())
    return;
  for (auto i : range(first, last)) {
    if (itvs[i].refine_as_integral()) {
      x.set_empty();
      return;
    }
  }
  if (keep_volume_info)
    x.maybe_update_volume_info();
  assert(x.check_inv());
}

template <bool KVI>
void
Box<KVI>::affine_image(Var var, Linear_Expr expr,
                       Integer inhomo, Integer den) {
  assert(space_dim() > 0);
  assert(var.space_dim() <= space_dim());
  assert(expr.space_dim() <= space_dim());
  assert(den != 0);

  if (is_empty())
    return;
  if (den < 0) {
    neg_assign(expr);
    neg_assign(inhomo);
    neg_assign(den);
  }
  assert(den > 0);

  auto vi = var.id();

  Index_Set pos, neg;
  auto get_sign_vars = [&expr, &pos, &neg]() {
    for (auto i : bwd_dim_range(expr) ) {
      auto s = sgn(expr.get(i));
      if (s == 0)
        continue;
      (s > 0) ?  pos.set(i) : neg.set(i);
    }
  };

  auto compute_new_bound =
    [this, &expr, &inhomo, &den](const Index_Set& lbs,
                                 const Index_Set& ubs) {
    Rational val(inhomo);
    for (auto i : lbs)
      add_mul_assign(val, Rational(expr[i]), itvs[i].lb);
    for (auto i : ubs)
      add_mul_assign(val, Rational(expr[i]), itvs[i].ub);
    val /= Rational(den);
    return val;
  };

  auto is_inf_lb = [this](dim_type i) { return itvs[i].inf_lb(); };
  auto is_inf_ub = [this](dim_type i) { return itvs[i].inf_ub(); };
  bool becomes_inf_lb = false;
  bool becomes_inf_ub = false;

  get_sign_vars();

  if (any_of(pos, is_inf_lb) || any_of(neg, is_inf_ub)) {
    becomes_inf_lb = true;
    unconstrain_lb(vi);
  }
  if (any_of(pos, is_inf_ub) || any_of(neg, is_inf_lb)) {
    becomes_inf_ub = true;
    unconstrain_ub(vi);
  }

  if (!becomes_inf_lb)
    itvs[vi].set_lb(compute_new_bound(pos, neg));

  if (!becomes_inf_ub)
    itvs[vi].set_ub(compute_new_bound(neg, pos));

  maybe_update_volume_info();
  assert(check_inv());
}

} // namespace pplite

#endif // !defined(pplite_BBox_impl_hh)
