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

#ifndef pplite_U_Poly_hh
#define pplite_U_Poly_hh 1

#include "globals.hh"
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
class U_Wrap {
public:
  using Kernel = PH;
  using Shell = U_Wrap<Kernel>;

  explicit
  U_Wrap(dim_type d = 0, Spec_Elem s = Spec_Elem::UNIVERSE,
         Topol t = get_default_topology())
    : info(d, not_a_dim()), kernel(0, s, t) { assert(check_inv()); }
  U_Wrap(dim_type d, Topol t, Spec_Elem s = Spec_Elem::UNIVERSE)
    : U_Wrap(d, s, t) {}
  U_Wrap(Spec_Elem s, dim_type d, Topol t = get_default_topology())
    : U_Wrap(d, s, t) {}
  U_Wrap(Spec_Elem s, Topol t, dim_type d)
    : U_Wrap(d, s, t) {}
  U_Wrap(Topol t, dim_type d, Spec_Elem s = Spec_Elem::UNIVERSE)
    : U_Wrap(d, s, t) {}
  U_Wrap(Topol t, Spec_Elem s, dim_type d)
    : U_Wrap(d, s, t) {}

  explicit U_Wrap(const Kernel& ph)
    : U_Wrap(ph.space_dim(), ph.topology(),
             ph.is_empty() ? Spec_Elem::EMPTY : Spec_Elem::UNIVERSE) {
    if (is_empty())
      return;
    const auto& cs = ph.cons();
    add_shell_cons(cs.begin(), cs.end());
  }

  // Implicit conversion to Kernel.
  operator Kernel() const {
    Kernel res = Kernel(space_dim(), topology(),
                        is_empty() ? Spec_Elem::EMPTY : Spec_Elem::UNIVERSE);
    if (!is_empty())
      res.add_cons(copy_cons());
    return res;
  }

  U_Wrap(const U_Wrap& y) = default;
  U_Wrap& operator=(const U_Wrap& y) = default;
  U_Wrap(U_Wrap&& y) noexcept = default;
  U_Wrap& operator=(U_Wrap&& y) noexcept = default;
  ~U_Wrap() = default;

  /* Predicates */
  bool check_inv() const;
  bool is_necessarily_closed() const { return kernel.is_necessarily_closed(); }
  bool is_empty() const { return kernel.is_empty(); }
  bool is_minimized() const { return kernel.is_minimized(); }
  bool is_universe() const {
    return !is_empty() && kernel.space_dim() == 0;
  }
  bool is_topologically_closed() const {
    return kernel.is_topologically_closed();
  }
  bool is_bounded() const {
    if (is_empty())
      return true;
    return space_dim() == kernel.space_dim() && kernel.is_bounded();
  }
  bool is_bounded_expr(bool from_below, const Linear_Expr& expr) const;
  bool constrains(Var var) const {
    return is_empty() || is_kernel_dim(var.id());
  }
  bool equals(const U_Wrap& y) const;
  bool contains(const U_Wrap& y) const;
  bool strictly_contains(const U_Wrap& y) const;
  bool is_disjoint_from(const U_Wrap& y) const;

  BBox get_bounding_box() const;
  bool boxed_contains(const U_Wrap& y) const;

  /* Queries */
  Topol topology() const { return kernel.topology(); }
  dim_type space_dim() const { return info.size(); }
  dim_type affine_dim() const {
    if (is_empty()) return 0;
    return space_dim() - kernel.space_dim() + kernel.affine_dim();
  }
  Poly_Con_Rel relation_with(const Con& c) const;
  Poly_Gen_Rel relation_with(const Gen& g) const;
  bool min(const Affine_Expr& ae, Rational& value,
           bool* included_ptr = nullptr, Gen* g_ptr = nullptr) const;
  bool max(const Affine_Expr& ae, Rational& value,
           bool* included_ptr = nullptr, Gen* g_ptr = nullptr) const {
    bool res = min(-ae, value, included_ptr, g_ptr);
    neg_assign(value);
    return res;
  }
  Itv get_bounds(Var var) const {
    if (is_empty())
      return Itv(Spec_Elem::EMPTY);
    auto sd = var.id();
    if (is_unconstrained_dim(sd))
      return Itv();
    auto k_var = Var(kernel_dim(sd));
    return kernel.get_bounds(k_var);
  }
  Itv get_bounds(const Affine_Expr& ae) const {
    if (is_empty())
      return Itv(Spec_Elem::EMPTY);
    for (auto sd : ae.expr.non_zeroes()) {
      if (is_unconstrained_dim(sd))
        return Itv();
    }
    auto k_ae = Affine_Expr(shell_2_kernel_expr(ae.expr), ae.inhomo);
    return kernel.get_bounds(k_ae);
  }
  Itv get_bounds(const Itv_Expr& ie) const {
    if (is_empty())
      return Itv(Spec_Elem::EMPTY);

    const auto& ie_vars = ie.first;
    const auto& ie_itvs = ie.second;
    assert(num_rows(ie_vars) == num_rows(ie_itvs));

    if (ie_vars.empty())
      return Itv::zero();

    // Check for trivial unboundedness.
    for (auto i : index_range(ie_vars)) {
      if (is_unconstrained_dim(ie_vars[i].id())
          && not ie_itvs[i].is_zero())
        return Itv();
    }

    Itv_Expr k_ie;
    Vars& k_vars = k_ie.first;
    Itvs& k_itvs = k_ie.second;
    for (auto i : index_range(ie_vars)) {
      auto ie_dim = ie_vars[i].id();
      if (is_kernel_dim(ie_dim)) {
        auto k_dim = kernel_dim(ie_dim);
        k_vars.push_back(Var(k_dim));
        k_itvs.push_back(ie_itvs[i]);
      }
    }
    return kernel.get_bounds(k_ie);
  }

  Index_Set get_unconstrained() const {
    Index_Set res;
    if (is_empty())
      return res;
    check_kernel_for_unconstrained();
    for (auto i : bwd_index_range(info)) {
      if (is_unconstrained_dim(i))
        res.set(i);
    }
    return res;
  }

  Cons copy_cons() const;
  Gens copy_gens() const;
  Cons copy_normalized_cons() const { return kernel.copy_normalized_cons(); }
  Gens_Info gens_info() const {
    auto info = kernel.gens_info();
    std::get<0>(info) += (space_dim() - kernel.space_dim());
    return info;
  }
  dim_type num_min_cons() const { return kernel.num_min_cons(); }
  dim_type num_min_gens() const {
    if (is_empty())
      return 0;
    return kernel.num_min_gens() + space_dim() - kernel.space_dim();
  }

  // A proxy type for constraints: maps kernel cs into shell cs.
  struct Cons_Proxy {
    using K_Proxy = typename Kernel::Cons_Proxy;
    K_Proxy k_proxy;
    Dims inverse;
    using value_type = Con;
    using cache_type = std::vector<std::unique_ptr<value_type>>;
    mutable cache_type shell_cs;

    explicit Cons_Proxy(const Shell* shell_ptr, bool normalized = false)
      : k_proxy(normalized ?
                shell_ptr->kernel.normalized_cons()
                : shell_ptr->kernel.cons()),
        inverse(shell_ptr->get_inverse_info()),
        shell_cs(k_proxy.end_pos()) {}

    // Movable.
    Cons_Proxy(Cons_Proxy&&) noexcept = default;
    Cons_Proxy& operator=(Cons_Proxy&&) noexcept = default;
    ~Cons_Proxy() = default;
    // Non-copiable.
    Cons_Proxy() = delete;
    Cons_Proxy(const Cons_Proxy&) = delete;
    Cons_Proxy& operator=(const Cons_Proxy&) = delete;

    bool is_skippable(dim_type pos) const { return k_proxy.is_skippable(pos); }
    dim_type end_pos() const { return num_rows(shell_cs); }
    dim_type size() const { return k_proxy.size(); }

    const value_type* get_value_ptr(dim_type pos) const {
      if (shell_cs[pos] == nullptr) {
        const auto& k_c = *(k_proxy.get_value_ptr(pos));
        auto s_c = Shell::kernel_2_shell_con(k_c, inverse);
        shell_cs[pos].reset(new value_type(std::move(s_c)));
      }
      return shell_cs[pos].get();
    }

    using const_iterator = Proxy_Iter<Cons_Proxy>;
    const_iterator cbegin() const { return const_iterator(this, false); }
    const_iterator cend() const { return const_iterator(this, true); }
    const_iterator begin() const { return cbegin(); }
    const_iterator end() const { return cend(); }
  }; // Cons_Proxy

  Cons_Proxy cons() const { return Cons_Proxy(this); }
  Cons_Proxy normalized_cons() const { return Cons_Proxy(this, true); }

  // A proxy type for generators: maps kernel gs into shell gs.
  struct Gens_Proxy {
    using K_Proxy = typename Kernel::Gens_Proxy;
    K_Proxy k_proxy;
    dim_type k_offset;
    Dims inverse;
    using value_type = Gen;
    using cache_type = std::vector<std::unique_ptr<value_type>>;
    mutable cache_type shell_gs;

    explicit Gens_Proxy(const Shell* shell_ptr)
      : k_proxy(shell_ptr->kernel.gens()),
        k_offset(0) {
      if (shell_ptr->is_empty())
        return;
      const auto sd = shell_ptr->space_dim();
      for (auto i : range(sd))
        if (shell_ptr->is_unconstrained_dim(i))
          ++k_offset;
      inverse = shell_ptr->get_inverse_info();
      shell_gs.reserve(k_offset + k_proxy.size());
      for (auto i : range(sd))
        if (shell_ptr->is_unconstrained_dim(i))
          shell_gs.emplace_back(new Gen(line(Var(i))));
      shell_gs.resize(k_offset + k_proxy.size());
    }

    // Movable.
    Gens_Proxy(Gens_Proxy&&) noexcept = default;
    Gens_Proxy& operator=(Gens_Proxy&&) noexcept = default;
    ~Gens_Proxy() = default;
    // Non-copiable.
    Gens_Proxy() = delete;
    Gens_Proxy(const Gens_Proxy&) = delete;
    Gens_Proxy& operator=(const Gens_Proxy&) = delete;

    // Nothing is skippable, hence size() == end_pos().
    bool is_skippable(dim_type) const { return false; }
    dim_type end_pos() const { return num_rows(shell_gs); }
    dim_type size() const { return num_rows(shell_gs); }

    const value_type* get_value_ptr(dim_type pos) const {
      if (shell_gs[pos] == nullptr) {
        assert(pos >= k_offset);
        auto k_pos = pos - k_offset;
        const auto& k_g = *(k_proxy.get_value_ptr(k_pos));
        auto s_g = Shell::kernel_2_shell_gen(k_g, inverse);
        shell_gs[pos].reset(new value_type(std::move(s_g)));
      }
      return shell_gs[pos].get();
    }

    using const_iterator = Proxy_Iter<Gens_Proxy>;
    const_iterator cbegin() const { return const_iterator(this, false); }
    const_iterator cend() const { return const_iterator(this, true); }
    const_iterator begin() const { return cbegin(); }
    const_iterator end() const { return cend(); }
  }; // Gens_Proxy

  Gens_Proxy gens() const { return Gens_Proxy(this); }

  void collapse(dim_type) { /* nothing to do */ }
  dim_type num_disjuncts() const { return is_empty() ? 0 : 1; }
  Cons_Proxy disjunct_cons(dim_type n) const {
    (void) n;
    assert(n == 0);
    return cons();
  }
  bool geom_covers(const U_Wrap& y) const { return contains(y); }

  /* Modifiers */

  void m_swap(U_Wrap& y) noexcept {
    using std::swap;
    swap(info, y.info);
    kernel.m_swap(y.kernel);
  }

  void set_empty() {
    info.assign(info.size(), not_a_dim());
    set_kernel(Spec_Elem::EMPTY);
  }
  void set_universe() {
    info.assign(info.size(), not_a_dim());
    set_kernel(Spec_Elem::UNIVERSE);
  }
  void set_topology(Topol t) { kernel.set_topology(t); }

  void add_con(Con c) { add_shell_cons(&c, &c + 1); }
  // Note: no point in passing move_iterators.
  void add_cons(Cons cs) { add_shell_cons(cs.begin(), cs.end()); }
  template <typename Iter>
  void add_cons(Iter first, Iter last) { add_shell_cons(first, last); }

  void add_gen(Gen g) { add_shell_gens(&g, &g + 1); }
  // Note: no point in passing move_iterators.
  void add_gens(Gens gs) { add_shell_gens(gs.begin(), gs.end()); }
  template <typename Iter>
  void add_gens(Iter first, Iter last) { add_shell_gens(first, last); }

  void topological_closure_assign() {
    kernel.topological_closure_assign();
  }

  template <typename Iter>
  void unconstrain(Iter first, Iter last) {
    if (first == last)
      return;
    if (is_empty())
      return;
    unconstrain_kernel(first, last);
    check_kernel_for_unconstrained();
    assert(check_inv());
  }
  void unconstrain(Var var) {
    auto d = var.id();
    unconstrain(&d, &d + 1);
  }
  void unconstrain(const Index_Set& vars) {
    unconstrain(vars.begin(), vars.end());
  }

  void intersection_assign(const U_Wrap& y);
  void join_assign(const U_Wrap& y) { poly_hull_assign(y); }
  void poly_difference_assign(const U_Wrap& y) {
    auto& x = *this;
    x = detail::poly_difference(x, y);
  }

private:
  template <typename Func, typename ...OptBool>
  void hull_assign(const U_Wrap& y, Func k_huller, OptBool&&... boxed) {
    auto& x = *this;
    if (y.is_empty() || x.is_universe())
      return;
    if (x.is_empty() || y.is_universe()) {
      x = y;
      return;
    }
    assert(x.kernel.space_dim() > 0 && y.kernel.space_dim() > 0);

    auto p = missing_kernel_dims(x, y);
    const auto& y_uncon = p.first;
    const auto& x_uncon = p.second;

    if (y_uncon.size() == y.kernel.space_dim()) {
      // y would become universe, hence also x.
      x.set_universe();
      return;
    }

    if (!x_uncon.empty()) {
      // Refrain from directly calling `unconstrain'
      // (since it may lead to kernels having different dimensions).
      x.unconstrain_kernel(x_uncon);
      // Check if we unconstrained all kernel dims.
      if (x.kernel.is_universe()) {
        x.set_universe();
        return;
      }
    }

    if (y_uncon.empty()) {
      x.sync_kernel_dims(y);
      k_huller(x.kernel, y.kernel, boxed...);
      x.check_kernel_for_unconstrained();
      return;
    }

    // Work on a copy of y.
    auto y_copy(y);
    // Refrain from directly calling `unconstrain' (ditto).
    y_copy.unconstrain_kernel(y_uncon);
    if (y_copy.kernel.is_universe()) {
      x.set_universe();
      return;
    }

    x.sync_kernel_dims(y_copy);
    k_huller(x.kernel, y_copy.kernel, boxed...);
    x.check_kernel_for_unconstrained();
    assert(check_inv());
  }

public:
  void poly_hull_assign(const U_Wrap& y) {
    hull_assign(y, std::mem_fn(&Kernel::poly_hull_assign));
  }
  void con_hull_assign(const U_Wrap& y, bool boxed = false) {
    hull_assign(y, std::mem_fn(&Kernel::con_hull_assign), boxed);
  }

  void affine_image(Var var, const Linear_Expr& expr,
                    const Integer& inhomo = Integer::zero(),
                    const Integer& den = Integer::one()) {
    affine_image_aux(false, var, expr, inhomo, den);
  }
  void affine_preimage(Var var, const Linear_Expr& expr,
                       const Integer& inhomo = Integer::zero(),
                       const Integer& den = Integer::one()) {
    affine_image_aux(true, var, expr, inhomo, den);
  }
  void
  parallel_affine_image(const Vars& vars,
                        const Linear_Exprs& exprs,
                        const Integers& inhomos, const Integers& dens);

  void widening_assign(const U_Wrap& y, const Cons* upto_ptr,
                       Widen_Impl w_impl, Widen_Spec w_spec);
  void widening_assign(const U_Wrap& y,
                       Widen_Impl w_impl = get_widen_impl(),
                       Widen_Spec w_spec = get_widen_spec()) {
    widening_assign(y, nullptr, w_impl, w_spec);
  }
  void widening_assign(const U_Wrap& y, const Cons& upto_cons,
                       Widen_Impl w_impl = get_widen_impl(),
                       Widen_Spec w_spec = get_widen_spec()) {
    widening_assign(y, &upto_cons, w_impl, w_spec);
  }

  void time_elapse_assign(const U_Wrap& y);

  /* Change of space dim */
  void add_space_dims(dim_type m, bool project = false);
  void map_space_dims(const Dims& pfunc);
  void remove_higher_space_dims(dim_type new_dim) {
    assert(new_dim <= space_dim());
    Dims vars(space_dim() - new_dim);
    std::iota(vars.begin(), vars.end(), new_dim);
    remove_space_dims(vars.begin(), vars.end());
  }
  template <typename Iter>
  void remove_space_dims(Iter first, Iter last) {
    unconstrain(first, last);
    erase_using_sorted_indices(info, first, last);
    assert(check_inv());
  }
  void remove_space_dims(const Index_Set& vars) {
    remove_space_dims(vars.begin(), vars.end());
  }
  void remove_space_dim(Var var) {
    auto dim = var.id();
    remove_space_dims(&dim, &dim + 1);
  }
  void concatenate_assign(const U_Wrap& y);
  void expand_space_dim(Var var, dim_type m);
  void fold_space_dims(const Index_Set& vars, Var dest);

private:
  U_Wrap split_aux(const Con& c, Topol t, bool integral);
public:
  U_Wrap split(const Con& c, Topol t) { return split_aux(c, t, false); }
  U_Wrap split(const Con& c) { return split(c, topology()); }
  U_Wrap integral_split(const Con& c) {
    return split_aux(c, topology(), true);
  }

  // Semantically const, but may affect syntactic representation.
  void minimize() const { kernel.minimize(); }

  size_t hash() const;

  /* Input-output */
  void print(std::ostream& s) const;
  void print() const { print(std::cout); }
  bool ascii_load(std::istream& s);
  void ascii_dump(std::ostream& s) const;

private:
  Dims info;
  Kernel kernel;

  Dims get_inverse_info() const;
  Cons kernel_2_shell_cons() const;
  Gens kernel_2_shell_gens() const;

  static Con kernel_2_shell_con(const Con& k_c, const Dims& inverse);
  static Gen kernel_2_shell_gen(const Gen& k_g, const Dims& inverse);

  bool is_kernel_dim(dim_type shell_dim) const {
    return info[shell_dim] != not_a_dim();
  }
  bool is_unconstrained_dim(dim_type shell_dim) const {
    return !is_kernel_dim(shell_dim);
  }

  dim_type& kernel_dim(dim_type shell_dim) {
    assert(is_kernel_dim(shell_dim));
    return info[shell_dim];
  }
  dim_type kernel_dim(dim_type shell_dim) const {
    assert(is_kernel_dim(shell_dim));
    return info[shell_dim];
  }

  bool all_kernel_dims(const Linear_Expr& expr) const {
    for (auto i : bwd_dim_range(expr)) {
      if (is_unconstrained_dim(i) && !expr.get(i).is_zero())
        return false;
    }
    return true;
  }

  static std::pair<Vars_Set, Vars_Set>
  missing_kernel_dims(const U_Wrap& x, const U_Wrap& y);

  Linear_Expr shell_2_kernel_expr(const Linear_Expr& s_expr) const;

  Con shell_2_kernel_con(const Con& s_c) const;
  Gen shell_2_kernel_gen(const Gen& s_g) const;

  void decr_info_helper(const Index_Set& k_dims) {
    dim_type found = 0;
    for (auto k_dim : k_dims) {
      for (auto& kd : info) {
        if (kd != not_a_dim() && kd > k_dim - found)
         --kd;
      }
      ++found;
    }
  }

  void set_kernel(Spec_Elem kind) {
    kernel = Kernel(0, kind, kernel.topology());
  }

  template <typename Iter>
  void add_shell_cons(Iter first, Iter last) {
    if (first == last)
      return;
    // First pass: update info with dims that will enter the kernel.
    const dim_type old_k_dim = kernel.space_dim();
    dim_type new_k_dims = 0;
    for (Iter i = first; i != last; ++i) {
      const auto& c = *i;
      for (dim_type j : dim_range(c)) {
        if (is_kernel_dim(j) || c.coeff(Var(j)).is_zero())
          continue;
        info[j] = old_k_dim + new_k_dims;
        ++new_k_dims;
      }
    }
    kernel.add_space_dims(new_k_dims);

    // Second pass: translate shell cons into kernel cons.
    Cons kernel_cs;
    for (Iter i = first; i != last; ++i) {
      const auto& s_c = *i;
      Con k_c = shell_2_kernel_con(s_c);
      if (!k_c.linear_expr().is_zero())
        kernel_cs.push_back(std::move(k_c));
      else if (k_c.is_inconsistent()) {
        set_empty();
        return;
      }
    }
    kernel.add_cons(std::move(kernel_cs));
    kernel.minimize();
    if (kernel.is_empty())
      set_empty();
  }

  template <typename Iter>
  void add_shell_gens(const Iter first, const Iter last) {
    if (first == last)
      return;
    if (is_universe())
      return;
    if (is_empty()) {
      std::iota(info.begin(), info.end(), 0);
      kernel.add_space_dims(space_dim());
      kernel.add_gens(first, last);
      check_kernel_for_unconstrained();
      assert(check_inv());
      return;
    }

    bool nothing_todo = true;
    for (Iter i = first; nothing_todo && i != last; ++i) {
      const auto& sg = *i;
      if (sg.is_point() || sg.is_closure_point()) {
        nothing_todo = false;
        break;
      }
      for (dim_type j : dim_range(sg)) {
        if (is_kernel_dim(j) && sg.coeff(Var(j)) != 0) {
          nothing_todo = false;
          break;
        }
      }
    }
    if (nothing_todo)
      return;

    Gens kernel_gs;
    for (Iter i = first; i != last; ++i) {
      const auto& s_g = *i;
      Linear_Expr k_expr;
      for (dim_type j : dim_range(s_g)) {
        if (is_kernel_dim(j) && s_g.coeff(Var(j)) != 0)
          k_expr.set(Var(kernel_dim(j)), s_g.coeff(Var(j)));
      }
      if ((s_g.is_ray() || s_g.is_line()) && k_expr.is_zero())
        continue;
      Gen k_g(s_g.type(), std::move(k_expr), s_g.impl().inhomo);
      kernel_gs.push_back(std::move(k_g));
    }
    kernel.add_gens(std::move(kernel_gs));
    kernel.minimize();
    check_kernel_for_unconstrained();
    assert(check_inv());
  }

  template <typename Iter>
  void unconstrain_kernel(Iter first, Iter last) {
    assert(first != last);
    assert(!is_empty());
    Index_Set k_dims;
    for ( ; first != last; ++first) {
      dim_type sd = *first;
      assert(sd < space_dim());
      if (is_kernel_dim(sd)) {
        k_dims.set(kernel_dim(sd));
        info[sd] = not_a_dim();
      }
    }
    // Remove the kernel variables (all at once).
    kernel.remove_space_dims(k_dims);
    // Fix info to match the removed kernel variables (one at a time).
    decr_info_helper(k_dims);
    // Note: NOT checking if other kernel dims became unconstrained.
    // This is up to the caller (class invariant may be broken here).
  }
  void unconstrain_kernel(const Index_Set& vars) {
    unconstrain_kernel(vars.begin(), vars.end());
  }

  void sync_kernel_dims(const U_Wrap& y) const;
  void sync_kernel_dims(const Dims& sync_info) const;
  void check_kernel_for_unconstrained() const;
  template <typename Fun>
  void sync_and_modify(const U_Wrap& y, Fun fun);

  void affine_image_aux(bool preimage, Var var, const Linear_Expr& expr,
                        const Integer& inhomo, const Integer& den);

}; // class U_Wrap

namespace IO_Operators {

template <typename PH>
inline std::ostream&
operator<<(std::ostream& s, const U_Wrap<PH>& ph) {
  ph.print(s);
  return s;
}

} // namespace IO_Operators

template <typename PH>
inline void
swap(U_Wrap<PH>& x, U_Wrap<PH>& y) noexcept { x.m_swap(y); }

template <typename PH>
inline bool
operator==(const U_Wrap<PH>& x, const U_Wrap<PH>& y) { return x.equals(y); }
template <typename PH>
inline bool
operator!=(const U_Wrap<PH>& x, const U_Wrap<PH>& y) { return !(x == y); }

// Explicit template instantiation (declaration).
extern template class U_Wrap<Poly>;
extern template class U_Wrap<F_Poly>;

using U_Poly = U_Wrap<Poly>;
using UF_Poly = U_Wrap<F_Poly>;

NOTHROW_MOVES(U_Poly);
NOTHROW_MOVES(UF_Poly);

} // namespace pplite

#endif // !defined(pplite_U_Poly_hh)
