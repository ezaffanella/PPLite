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

#include "pplite-config.h"
#include "U_Poly.hh"
#include "F_Poly.hh"
#include "Poly_widen.hh"

#include <set>
#include <iostream>

namespace pplite {

template <typename PH>
bool
U_Wrap<PH>::check_inv() const {
  std::string reason;
#ifdef NDEBUG
  auto maybe_dump = []() {};
#else // In debugging mode, be noisy.
  auto maybe_dump = [this, &reason]() {
    std::cerr << reason << std::endl;
    this->ascii_dump(std::cerr);
    std::cerr << reason << std::endl;
  };
#endif
  // Checks on info.
  dim_type num_kdims = 0;
  std::set<dim_type> found_kernel_dim;
  for (auto sd : bwd_dim_range(*this)) {
    if (is_kernel_dim(sd)) {
      ++num_kdims;
      const auto k_dim = kernel_dim(sd);
      if (k_dim >= kernel.space_dim()) {
        reason = "U_Wrap broken: KERNEL dimension not in kernel";
        maybe_dump();
        return false;
      }
      if (found_kernel_dim.find(k_dim) != found_kernel_dim.end()) {
        reason = "U_Wrap broken: info is not injective";
        maybe_dump();
        return false;
      }
      found_kernel_dim.insert(k_dim);
    }
  }
  if (num_kdims != kernel.space_dim()) {
    reason = "U_Wrap broken: kernel space dimension does not match "
      "the number of KERNEL dimensions";
    maybe_dump();
    return false;
  }
  // Checks on kernel.
  if (is_empty() && kernel.space_dim() > 0) {
    reason = "U_Wrap kernel broken: empty polyhedron has KERNEL dimensions";
    maybe_dump();
    return false;
  }
  for (auto i : bwd_dim_range(kernel) ) {
    if (!kernel.constrains(Var(i))) {
      reason = "U_Wrap kernel broken: kernel has unconstrained dimension";
      maybe_dump();
      return false;
    }
  }
  // All checks passed.
  return true;
}

template <typename PH>
Dims
U_Wrap<PH>::get_inverse_info() const {
  // Build inverse map (kdim -> sdim).
  const auto k_dim = kernel.space_dim();
  Dims inverse(k_dim);
  for (auto i : bwd_dim_range(*this)) {
    if (is_kernel_dim(i))
      inverse[kernel_dim(i)] = i;
  }
  return inverse;
}

template <typename PH>
Linear_Expr
U_Wrap<PH>::shell_2_kernel_expr(const Linear_Expr& s_expr) const {
  Linear_Expr k_expr;
  for (auto i : dim_range(s_expr)) {
    if (is_kernel_dim(i))
      k_expr.set(kernel_dim(i), s_expr.get(i));
  }
  return k_expr;
}

template <typename PH>
inline Con
U_Wrap<PH>::shell_2_kernel_con(const Con& s_c) const {
  Linear_Expr k_expr = shell_2_kernel_expr(s_c.linear_expr());
  return Con(std::move(k_expr), s_c.inhomo_term(), s_c.type());
}

template <typename PH>
inline Gen
U_Wrap<PH>::shell_2_kernel_gen(const Gen& s_g) const {
  Linear_Expr k_expr = shell_2_kernel_expr(s_g.linear_expr());
  assert(not (s_g.is_line_or_ray() && k_expr.is_zero()));
  return Gen(s_g.type(), std::move(k_expr), s_g.impl().inhomo);
}

template <typename PH>
void
U_Wrap<PH>::check_kernel_for_unconstrained() const {
  // Note: all changes preserve semantics.
  auto& x = const_cast<U_Wrap&>(*this);
  Vars_Set k_remove;
  for (auto sd : bwd_dim_range(x)) {
    if (x.is_kernel_dim(sd)) {
      const Var k_var(x.kernel_dim(sd));
      if (!x.kernel.constrains(k_var)) {
        k_remove.insert(k_var);
        x.info[sd] = not_a_dim();
      }
    }
  }
  x.kernel.remove_space_dims(k_remove);
  x.kernel.minimize();
  x.decr_info_helper(k_remove);
}

template <typename PH>
Con
U_Wrap<PH>::kernel_2_shell_con(const Con& k_c,
                               const Dims& inverse) {
  const auto& k_expr = k_c.linear_expr();
  Linear_Expr s_expr;
  for (auto j : bwd_dim_range(k_expr)) {
    if (!k_expr[j].is_zero())
      s_expr.set(Var(inverse[j]), k_expr[j]);
  }
  return Con(std::move(s_expr), k_c.inhomo_term(), k_c.type());
}

template <typename PH>
Cons
U_Wrap<PH>::copy_cons() const {
  Cons result;
  Dims inverse = get_inverse_info();
  for (const auto& k_c : kernel.cons())
    result.push_back(kernel_2_shell_con(k_c, inverse));
  return result;
}

template <typename PH>
Gen
U_Wrap<PH>::kernel_2_shell_gen(const Gen& k_g,
                               const Dims& inverse) {
  Linear_Expr s_expr;
  const auto& k_expr = k_g.linear_expr();
  for (auto i : dim_range(k_expr)) {
    if (!k_expr[i].is_zero())
      s_expr.set(Var(inverse[i]), k_expr[i]);
  }
  return Gen(k_g.type(), std::move(s_expr), k_g.impl().inhomo);
}

template <typename PH>
Gens
U_Wrap<PH>::copy_gens() const {
  Gens result;
  if (is_empty())
    return result;
  for (auto i : dim_range(*this)) {
    if (!is_kernel_dim(i))
      result.push_back(line(Var(i)));
  }
  Dims inverse = get_inverse_info();
  for (const auto& k_g : kernel.gens())
    result.push_back(kernel_2_shell_gen(k_g, inverse));
  return result;
}

template <typename PH>
void
U_Wrap<PH>::sync_kernel_dims(const Dims& sync_info) const {
  // Note: all changes preserve semantics.
  auto& x = const_cast<U_Wrap&>(*this);
  assert(x.space_dim() == num_rows(sync_info));
  // Permute the kernel dims of x to match the info in sync.
  Dims mapper(x.kernel.space_dim(), not_a_dim());
  for (auto sd : bwd_dim_range(x)) {
    if (x.is_kernel_dim(sd)) {
      auto new_dim = sync_info[sd];
      assert(new_dim != not_a_dim());
      mapper[x.kernel_dim(sd)] = new_dim;
      x.kernel_dim(sd) = new_dim;
    }
  }
  assert(x.info == sync_info);
  x.kernel.map_space_dims(mapper);
}

template <typename PH>
void
U_Wrap<PH>::sync_kernel_dims(const U_Wrap& y) const {
  assert(space_dim() == y.space_dim());
  assert(kernel.space_dim() == y.kernel.space_dim());
  sync_kernel_dims(y.info);
}

template <typename PH>
Poly_Con_Rel
U_Wrap<PH>::relation_with(const Con& c) const {
  assert(space_dim() >= c.space_dim());

  if (is_empty())
    return Poly_Con_Rel::saturates()
      && Poly_Con_Rel::is_included()
      && Poly_Con_Rel::is_disjoint();
  if (space_dim() == 0)
    return kernel.relation_with(c);
  check_kernel_for_unconstrained();
  if (!all_kernel_dims(c.linear_expr()))
    return Poly_Con_Rel::strictly_intersects();
  return kernel.relation_with(shell_2_kernel_con(c));
}

template <typename PH>
Poly_Gen_Rel
U_Wrap<PH>::relation_with(const Gen& g) const {
  assert(space_dim() >= g.space_dim());
  if (is_empty())
    return Poly_Gen_Rel::nothing();
  if (space_dim() == 0)
    return Poly_Gen_Rel::subsumes();
  if (g.is_line() || g.is_ray()) {
    // Quick check for subsumption.
    bool quick_subsumes = true;
    for (auto i : bwd_dim_range(g)) {
      if (is_kernel_dim(i) && g.coeff(Var(i)) != 0) {
        quick_subsumes = false;
        break;
      }
    }
    if (quick_subsumes)
      return Poly_Gen_Rel::subsumes();
  }
  return kernel.relation_with(shell_2_kernel_gen(g));
}

template <typename PH>
bool
U_Wrap<PH>::equals(const U_Wrap& y) const {
  const auto& x = *this;
  if (x.kernel.space_dim() != y.kernel.space_dim())
    return false;
  for (auto sd : bwd_dim_range(x)) {
    if (x.is_kernel_dim(sd) ^ y.is_kernel_dim(sd))
      return false;
  }
  x.sync_kernel_dims(y);
  return x.kernel.equals(y.kernel);
}

template <typename PH>
bool
U_Wrap<PH>::contains(const U_Wrap& y) const {
  const auto& x = *this;
  if (x.is_empty())
    return y.is_empty();
  if (y.is_empty() || x.is_universe())
    return true;

  auto p = missing_kernel_dims(x, y);
  const auto& x_missing = p.first;
  const auto& y_missing = p.second;

  if (!y_missing.empty())
    return false;
  if (x_missing.size() == y.kernel.space_dim())
    return false;
  if (x_missing.empty()) {
    x.sync_kernel_dims(y);
    return x.kernel.contains(y.kernel);
  }
  // Work on a copy of y.
  auto y_copy(y);
  y_copy.unconstrain_kernel(x_missing);
  assert(x.kernel.space_dim() == y_copy.kernel.space_dim());
  x.sync_kernel_dims(y_copy);
  return x.kernel.contains(y_copy.kernel);
}

template <typename PH>
bool
U_Wrap<PH>::strictly_contains(const U_Wrap& y) const {
  const auto& x = *this;
  if (x.is_empty())
    return false;
  if (y.is_empty())
    return true;
  if (x.is_universe() && !y.is_universe())
    return true;

  auto p = missing_kernel_dims(x, y);
  const auto& x_missing = p.first;
  const auto& y_missing = p.second;

  if (!y_missing.empty())
    return false;
  if (x_missing.size() == y.kernel.space_dim())
    return false;
  if (x_missing.empty()) {
    x.sync_kernel_dims(y);
    return x.kernel.strictly_contains(y.kernel);
  }
  // Work on a copy of y.
  auto y_copy(y);
  y_copy.unconstrain_kernel(x_missing);
  assert(x.kernel.space_dim() == y_copy.kernel.space_dim());
  x.sync_kernel_dims(y_copy);
  // Note: since we unconstrained some dims of y, we answer positively
  // even if kernel containment is non-strict.
  return x.kernel.contains(y_copy.kernel);
}

template <typename PH>
bool
U_Wrap<PH>::is_disjoint_from(const U_Wrap& y) const {
  const auto& x = *this;
  if (x.is_empty() || y.is_empty())
    return true;

  auto p = missing_kernel_dims(x, y);
  const auto& x_missing = p.first;
  const auto& y_missing = p.second;

  auto uncon_copy = [](const U_Wrap& ph, const Vars_Set& uncon_vs) {
    auto ph_copy(ph);
    ph_copy.unconstrain_kernel(uncon_vs);
    return ph_copy;
  };

  if (x_missing.empty() && y_missing.empty()) {
    x.sync_kernel_dims(y);
    return x.kernel.is_disjoint_from(y.kernel);
  }
  if (x_missing.empty()) {
    assert(!y_missing.empty());
    auto x_copy = uncon_copy(x, y_missing);
    x_copy.sync_kernel_dims(y);
    return x_copy.kernel.is_disjoint_from(y.kernel);
  }
  if (y_missing.empty()) {
    assert(!x_missing.empty());
    auto y_copy = uncon_copy(y, x_missing);
    y_copy.sync_kernel_dims(x);
    return x.kernel.is_disjoint_from(y_copy.kernel);
  }
  assert(!x_missing.empty());
  assert(!y_missing.empty());
  auto x_copy = uncon_copy(x, y_missing);
  auto y_copy = uncon_copy(y, x_missing);
  x_copy.sync_kernel_dims(y_copy);
  return x_copy.kernel.is_disjoint_from(y_copy.kernel);
}

template <typename PH>
BBox
U_Wrap<PH>::get_bounding_box() const {
  if (is_empty())
    return BBox(space_dim(), Spec_Elem::EMPTY);
  auto kernel_res = kernel.get_bounding_box();
  auto res = BBox(space_dim());
  auto& rays = res.volume.first;
  auto& pvol = res.volume.second;
  rays = 2 * (space_dim() - kernel.space_dim()) + kernel_res.volume.first;
  pvol = (rays == 0) ? kernel_res.volume.second : Rational::zero();
  // Move kernel bounds into shell bounds.
  for (auto sd : bwd_dim_range(*this)) {
    if (is_kernel_dim(sd)) {
      auto kd = kernel_dim(sd);
      res[sd] = std::move(kernel_res[kd]);
    }
  }
  return res;
}

template <typename PH>
bool
U_Wrap<PH>::boxed_contains(const U_Wrap& y) const {
  const auto& x = *this;
  if (x.space_dim() == 0 || x.is_empty() || y.is_empty())
    return true;
  if (x.is_universe())
    return true;
  // Limit check to x cons that are NOT closed interval cons.
  for (const auto& c : x.cons()) {
    if (!c.is_strict_inequality() && is_proper_interval_con(c))
      continue;
    if (!y.relation_with(c).implies(Poly_Con_Rel::is_included()))
      return false;
  }
  return true;
}

template <typename PH>
std::pair<Vars_Set, Vars_Set>
U_Wrap<PH>::missing_kernel_dims(const U_Wrap& x, const U_Wrap& y) {
  assert(x.space_dim() == y.space_dim());
  Vars_Set x_res;
  Vars_Set y_res;
  for (auto sd : bwd_dim_range(x)) {
    if (!x.is_kernel_dim(sd)) {
      if (y.is_kernel_dim(sd))
        x_res.set(sd);
      continue;
    }
    if (!y.is_kernel_dim(sd))
      y_res.set(sd);
  }
  return std::make_pair(x_res, y_res);
}

namespace detail {

template <typename PH>
void
kernelize(Dims& info, PH& kernel, const Vars_Set& vars) {
  const dim_type old_dim = kernel.space_dim();
  dim_type new_dim = 0;
  for (auto sd : vars) {
    info[sd] = old_dim + new_dim;
    ++new_dim;
  }
  kernel.add_space_dims(new_dim);
}

} // namespace detail

template <typename PH>
template <typename Fun>
void
U_Wrap<PH>::sync_and_modify(const U_Wrap<PH>& y, Fun fun) {
  auto& x = *this;
  assert(x.space_dim() == y.space_dim());

  auto p = missing_kernel_dims(x, y);
  const auto& x_missing = p.first;
  detail::kernelize(x.info, x.kernel, x_missing);
  const auto& y_missing = p.second;
  if (y_missing.empty()) {
    x.sync_kernel_dims(y);
    fun(x.kernel, y.kernel);
  } else {
    // Work on a copy of y.
    auto y_copy = y;
    detail::kernelize(y_copy.info, y_copy.kernel, y_missing);
    x.sync_kernel_dims(y_copy);
    fun(x.kernel, y_copy.kernel);
  }
}

template <typename PH>
void
U_Wrap<PH>::intersection_assign(const U_Wrap& y) {
  auto& x = *this;
  if (x.is_empty() || y.is_universe())
    return;
  if (y.is_empty() || x.is_universe()) {
    x = y;
    return;
  }
  // Both non-empty and non-universe (hence, non zero-dim).
  assert(x.kernel.space_dim() > 0 && y.kernel.space_dim() > 0);

  x.sync_and_modify(y, std::mem_fn(&Kernel::intersection_assign));
  if (x.kernel.is_empty())
    x.set_empty();

  assert(check_inv());
}

template <typename PH>
void
U_Wrap<PH>::time_elapse_assign(const U_Wrap& y) {
  auto& x = *this;
  if (x.is_empty() || y.is_empty()) {
    x.set_empty();
    return;
  }
  if (x.is_universe() || y.is_universe()) {
    x.set_universe();
    return;
  }

  Index_Set uncon;
  for (auto sd : bwd_dim_range(x)) {
    if (y.is_unconstrained_dim(sd) && x.is_kernel_dim(sd))
      uncon.set(sd);
  }
  x.unconstrain(uncon);
  if (x.is_universe())
    return;

  Gens y_s_gens;
  Dims y_inv = y.get_inverse_info();
  for (const auto& k_g : y.kernel.gens())
    y_s_gens.push_back(kernel_2_shell_gen(k_g, y_inv));

  Gens rays;
  detail::add_as_rays(std::move(y_s_gens), rays);
  x.add_gens(std::move(rays));
  assert(x.check_inv());
}

template <typename PH>
void
U_Wrap<PH>::widening_assign(const U_Wrap& y, const Cons* upto_ptr,
                            Widen_Impl w_impl, Widen_Spec w_spec) {
  auto& x = *this;
  if (detail::widening_preamble(x, y, w_spec))
    return;
  if (x.is_universe())
    return;

  // Not a special case.
  if (w_spec == Widen_Spec::SAFE) {
    // Apply trivial lifting of risky widening.
    x.poly_hull_assign(y);
    x.minimize();
  }

  Index_Set valid = detail::valid_upto_cons(x, upto_ptr);
  // Apply (risky) widening to kernel poly.
  const auto p = missing_kernel_dims(x, y);
  const auto& x_missing = p.first;
  assert(p.second.empty());
  detail::kernelize(x.info, x.kernel, x_missing);
  x.sync_kernel_dims(y);
  x.kernel.widening_assign(y.kernel, w_impl, Widen_Spec::RISKY);
  x.check_kernel_for_unconstrained();

  detail::add_valid_upto_cons(x, valid, upto_ptr);

  assert(x.check_inv());
}

template <typename PH>
void
U_Wrap<PH>::concatenate_assign(const U_Wrap& y) {
  auto& x = *this;
  const dim_type old_sdim = x.space_dim();
  x.info.insert(x.info.end(), y.info.begin(), y.info.end());
  if (x.is_empty() || y.is_empty()) {
    x.set_empty();
    return;
  }
  // Adjust concatenated info (if needed).
  const dim_type old_kdim = x.kernel.space_dim();
  if (old_kdim > 0 && y.kernel.space_dim() > 0) {
    for (auto i : range(old_sdim, x.space_dim())) {
      if (x.is_kernel_dim(i))
        x.kernel_dim(i) += old_kdim;
    }
  }
  // Concatenate (non-empty) kernels.
  x.kernel.concatenate_assign(y.kernel);
  assert(check_inv());
}

template <typename PH>
void
U_Wrap<PH>::add_space_dims(dim_type m, bool project) {
  if (m == 0)
    return;
  const dim_type new_dim = space_dim() + m;
  if (!project || is_empty()) {
    info.resize(new_dim, not_a_dim());
    assert(check_inv());
    return;
  }
  const dim_type kd = kernel.space_dim();
  info.reserve(new_dim);
  for (auto i : range(m))
    info.push_back(kd + i);
  kernel.add_space_dims(m, project);
  // FIXME: needed? kernel.minimize();
  assert(check_inv());
}

template <typename PH>
void
U_Wrap<PH>::map_space_dims(const Dims& pfunc) {
  if (space_dim() == 0)
    return;
  assert(space_dim() == num_rows(pfunc));
  const auto max_dim = *std::max_element(pfunc.begin(), pfunc.end());
  if (max_dim == not_a_dim()) {
    // Empty codomain.
    info.clear();
    if (kernel.is_empty())
      set_empty();
    else
      set_universe();
    return;
  }

  const auto old_s_dim = space_dim();
  const auto new_s_dim = max_dim + 1;
  // First, unconstrain those kernel variables that will be removed.
  Vars_Set uncon;
  for (auto i : bwd_range(old_s_dim)) {
    if (is_kernel_dim(i) && pfunc[i] == not_a_dim())
      uncon.set(i);
  }
  unconstrain(uncon);
  // Then, compute new shell-kernel mapping.
  Dims new_info(new_s_dim, not_a_dim());
  for (auto i : bwd_range(old_s_dim)) {
    const auto j = pfunc[i];
    if (j != not_a_dim())
      new_info[j] = info[i];
  }
  using std::swap;
  swap(info, new_info);
  assert(check_inv());
}

template <typename PH>
void
U_Wrap<PH>::expand_space_dim(Var v, dim_type m) {
  if (m == 0)
    return;
  assert(m > 0);
  const dim_type sd = v.id();
  const dim_type new_sdim = space_dim() + m;
  if (is_unconstrained_dim(sd) || is_empty()) {
    info.resize(new_sdim, not_a_dim());
    assert(check_inv());
    return;
  }
  const dim_type old_kdim = kernel.space_dim();
  info.reserve(new_sdim);
  for (auto i : range(m))
    info.push_back(old_kdim + i);
  Var kv(kernel_dim(sd));
  kernel.expand_space_dim(kv, m);
  assert(check_inv());
}

template <typename PH>
void
U_Wrap<PH>::fold_space_dims(const Index_Set& vars, Var dest) {
  if (is_unconstrained_dim(dest.id())) {
    remove_space_dims(vars);
    return;
  }
  if (any_of(vars, [this](dim_type i) { return is_unconstrained_dim(i); })) {
    // vars contains an unconstrained dimension.
    remove_space_dims(vars);
    unconstrain(dest);
    assert(check_inv());
    return;
  }
  // All dims (in vars and dest) are constrained.
  Var k_dest(kernel_dim(dest.id()));
  Index_Set k_vars;
  for (auto i : vars) {
    k_vars.set(kernel_dim(i));
    info[i] = not_a_dim();
  }
  erase_using_sorted_indices(info, vars);
  decr_info_helper(k_vars);
  kernel.fold_space_dims(k_vars, k_dest);
  check_kernel_for_unconstrained();
  assert(check_inv());
}

template <typename PH>
U_Wrap<PH>
U_Wrap<PH>::split_aux(const Con& c, Topol t, bool integral) {
  assert(not c.is_equality() || integral);
  if (is_empty())
    return *this;
  minimize();
  const dim_type old_k_dim = kernel.space_dim();
  dim_type new_k_dims = 0;
  for (auto j : dim_range(c)) {
    if (is_kernel_dim(j) || c.coeff(Var(j)).is_zero())
      continue;
    info[j] = old_k_dim + new_k_dims;
    ++new_k_dims;
  }
  kernel.add_space_dims(new_k_dims);
  Con k_c = shell_2_kernel_con(c);

  auto res = U_Wrap(0, Spec_Elem::EMPTY, topology());
  res.info = info;
  res.kernel = integral
    ? kernel.integral_split(k_c)
    : kernel.split(k_c, t);

  if (kernel.is_empty())
    set_empty();
  if (res.kernel.is_empty())
    res.set_empty();
  else if (integral && c.is_equality())
    res.check_kernel_for_unconstrained();
  assert(check_inv() && res.check_inv());
  return res;
}

template <typename PH>
void
U_Wrap<PH>::affine_image_aux(bool preimage,
                             Var var, const Linear_Expr& expr,
                             const Integer& inhomo, const Integer& den) {
  if (is_empty())
    return;
  if (is_unconstrained_dim(var.id()) && expr.get(var) != 0)
    // Invertible affine (pre-)image on an unconstrained variable: nop.
    return;

  // Update info with dimensions that will enter the kernel.
  const dim_type old_k_dim = kernel.space_dim();
  dim_type new_k_dims = 0;
  if (is_unconstrained_dim(var.id())) {
    info[var.id()] = old_k_dim + new_k_dims;
    ++new_k_dims;
  }
  for (auto j : dim_range(expr)) {
    if (is_unconstrained_dim(j) && expr.get(j) != 0) {
      assert(j != var.id());
      info[j] = old_k_dim + new_k_dims;
      ++new_k_dims;
    }
  }
  kernel.add_space_dims(new_k_dims);
  // Rewrite var and expr for kernel.
  Var k_var(kernel_dim(var.id()));
  Linear_Expr k_expr = shell_2_kernel_expr(expr);
  // Apply affine (pre-)image on kernel.
  if (preimage)
    kernel.affine_preimage(k_var, k_expr, inhomo, den);
  else
    kernel.affine_image(k_var, k_expr, inhomo, den);
  check_kernel_for_unconstrained();
  assert(check_inv());
}

template <typename PH>
void
U_Wrap<PH>::parallel_affine_image(const Vars& vars,
                                  const Linear_Exprs& exprs,
                                  const Integers& inhomos,
                                  const Integers& dens) {
  const auto nvars = num_rows(vars);
  if (nvars == 1) {
    affine_image(vars[0], exprs[0], inhomos[0], dens[0]);
    return;
  }
  if (is_empty())
    return;

  // Update info with dimensions that will enter the kernel.
  const dim_type old_k_dim = kernel.space_dim();
  dim_type new_k_dims = 0;
  for (auto i : range(nvars)) {
    const auto sd = vars[i].id();
    if (is_unconstrained_dim(sd)) {
      info[sd] = old_k_dim + new_k_dims;
      ++new_k_dims;
    }
    const auto& e_i = exprs[i];
    for (auto j : dim_range(e_i)) {
      if (is_kernel_dim(j) || e_i.get(j).is_zero())
        continue;
      info[j] = old_k_dim + new_k_dims;
      ++new_k_dims;
    }
  }
  kernel.add_space_dims(new_k_dims);
  // Rewrite vars and exprs for kernel.

  Vars k_vars(nvars, Var(0));
  Linear_Exprs k_exprs(nvars);
  for (auto i : range(nvars)) {
    const auto sd = vars[i].id();
    k_vars[i] = Var(kernel_dim(sd));
    k_exprs[i] = shell_2_kernel_expr(exprs[i]);
  }
  // Apply parallel affine image on kernel.
  kernel.parallel_affine_image(k_vars, k_exprs, inhomos, dens);
  check_kernel_for_unconstrained();
  assert(check_inv());
}

template <typename PH>
bool
U_Wrap<PH>::is_bounded_expr(bool from_below, const Linear_Expr& expr) const {
  if (space_dim() == 0 || is_empty())
    return true;
  if (!all_kernel_dims(expr))
    return false;
  return kernel.is_bounded_expr(from_below, shell_2_kernel_expr(expr));
}

template <typename PH>
bool
U_Wrap<PH>::min(const Affine_Expr& ae, Rational& value,
                bool* included_ptr, Gen* g_ptr) const {
  if (is_empty())
    return false;
  if (!all_kernel_dims(ae.expr))
    return false;
  auto k_ae = Affine_Expr(shell_2_kernel_expr(ae.expr), ae.inhomo);
  auto has_min = kernel.min(k_ae, value, included_ptr, g_ptr);
  if (not has_min)
    return false;
  Dims inverse = get_inverse_info();
  if (g_ptr) {
    Gen g = kernel_2_shell_gen(*g_ptr, inverse);
    *g_ptr = std::move(g);
  }
  return true;
}

template <typename PH>
size_t
U_Wrap<PH>::hash() const {
  auto is_normalized = [](const Dims& info) {
    dim_type norm = 0;
    for (auto kd : info) {
      if (kd != not_a_dim()) {
        if (kd == norm)
          ++norm;
        else
          return false;
      }
    }
    return true;
  };
  auto get_normalized = [](const Dims& info) {
    Dims res = info;
    dim_type norm = 0;
    for (auto& kd : res) {
      if (kd != not_a_dim()) {
        kd = norm;
        ++norm;
      }
    }
    return res;
  };

  if (!is_normalized(info)) {
    auto norm_info = get_normalized(info);
    sync_kernel_dims(norm_info);
    assert(is_normalized(info));
  }
  size_t res = hash_size(space_dim());
  for (auto k : info)
    hash_combine(res, k);
  hash_combine(res, kernel.hash());
  return res;
}

template <typename PH>
void
U_Wrap<PH>::print(std::ostream& s) const {
  if (is_empty()) {
    s << "false";
    return;
  }
  minimize();
  bool comma = false;
  for (const auto& c : cons()) {
    if (comma) s << ", ";
    using namespace IO_Operators;
    s << c;
    comma = true;
  }
}

template <typename PH>
bool
U_Wrap<PH>::ascii_load(std::istream& s) {
  std::string str;
  dim_type s_dim;
  if (not (ascii_load_string(s, "space_dim") && (s >> s_dim) && (s_dim >= 0)))
    return false;
  if (!(s >> str) || (str != "empty" && str != "not_empty"))
    return false;
  bool flag_empty = (str == "empty");
  info.resize(s_dim, not_a_dim());
  for (auto i : range(s_dim)) {
    dim_type j;
    if (not (ascii_load_string(s, "dim")))
      return false;
    if (!(s >> j) || j != i)
      return false;
    if (not (ascii_load_string(s, "->")))
      return false;
    if (!(s >> str))
      return false;
    if (str != "not_a_dim")
      info[i] = std::stol(str);
  }
  if (!kernel.ascii_load(s))
    return false;
  if (kernel.is_empty() != flag_empty)
    return false;
  // Check invariants.
  assert(check_inv());
  return true;
}

template <typename PH>
void
U_Wrap<PH>::ascii_dump(std::ostream& s) const {
  s << "space_dim " << space_dim();
  s << " " << (is_empty() ? "empty" : "not_empty") << "\n";
  for (auto i : index_range(info)) {
    s << "dim " << i << " -> ";
    if (is_kernel_dim(i))
      s << kernel_dim(i);
    else
      s << "not_a_dim";
    s << "\n";
  }
  s << "\n";
  kernel.ascii_dump(s);
  s << "\n";
}

// Explicit template instantiation (definition).
template class U_Wrap<Poly>;
template class U_Wrap<F_Poly>;

} // namespace pplite
