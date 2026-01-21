/* PPLite: a lightweight library for convex polyhedra derived from PPL.
   Copyright (C) 2001-2010 Roberto Bagnara <bagnara@cs.unipr.it>
   Copyright (C) 2010-2018 BUGSENG srl (http://bugseng.com)
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

#include "pplite-config.h"
#include "PolySet.hh"
#include "PolySet_templ.hh"
#include "Poly_widen.hh"

#include <algorithm>
#include <functional>
#include <iostream>
#include <list>
#include <map>
#include <string>

namespace pplite {

template <typename PH>
bool
PolySet<PH>::check_omega_reduced(const Seq& seq) {
  auto first = seq.begin();
  auto last = seq.end();
  for (auto i = first; i != last; ++i) {
    if (i->is_empty())
      return false;
    auto j = i;
    for (++j; j != last; ++j) {
      if (i->contains(*j) || j->contains(*i))
        return false;
    }
  }
  return true;
}

template <typename PH>
inline void
PolySet<PH>::linear_partition_aux(const Con& c, DS_Pair& p) {
  auto& ph = p.first;
  auto& seq = p.second;
  auto ph_neg = ph.split(c, Topol::NNC);
  if (!ph_neg.is_empty())
    seq.push_back(std::move(ph_neg));
}

template <typename PH>
typename PolySet<PH>::DS_Pair
PolySet<PH>::linear_partition(const Disj& p, const Disj& q) {
  auto res = DS_Pair(q, Seq());
  for (const auto& c : p.cons()) {
    if (c.is_equality()) {
      // FIXME: three way split here?
      Con c1(c.linear_expr(), c.inhomo_term(), Con::NONSTRICT_INEQUALITY);
      linear_partition_aux(c1, res);
      Con c2(-c.linear_expr(), -c.inhomo_term(), Con::NONSTRICT_INEQUALITY);
      linear_partition_aux(c2, res);
    } else
      linear_partition_aux(c, res);
  }
  return res;
}

template <typename PH>
bool
PolySet<PH>::check_containment(const Disj& d, const Seq& seq) {
  if (d.is_empty())
    return true;

  Seq tmp_seq { d };
  for (const auto& di : seq) {
    tmp_seq.remove_if([&di](const Disj& dj) { return di.contains(dj); });
    if (tmp_seq.empty())
      return true;
    Seq new_seq;
    for (auto j = tmp_seq.begin(); j != tmp_seq.end(); ) {
      if (j->is_disjoint_from(di))
        ++j;
      else {
        auto p = linear_partition(di, *j);
        new_seq.splice(new_seq.end(), p.second);
        j = tmp_seq.erase(j);
      }
    }
    tmp_seq.splice(tmp_seq.end(), new_seq);
    // FIXME: omega reduction here?
  }
  return false;
}

template <typename PH>
bool
PolySet<PH>::check_inv() const {
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
  for (const auto& d : seq()) {
    if (d.space_dim() != space_dim()) {
      reason = "PolySet sequence broken: disjunct has wrong space dim";
      maybe_dump();
      return false;
    }
    if (d.topology() != topology()) {
      reason = "PolySet sequence broken: disjunct has wrong topology";
      maybe_dump();
      return false;
    }
    if (!d.check_inv()) {
      reason = "PolySet sequence broken: disjunct is broken";
      return false;
    }
  }
  if (is_reduced() && !check_omega_reduced(seq())) {
    reason = "PolySet broken: claims to be reduced, but it is not";
    maybe_dump();
    return false;
  }
  // All checks passed.
  return true;
}

template <typename PH>
bool
PolySet<PH>::equals(const PolySet& y) const {
  const auto& x = *this;
  x.omega_reduce();
  y.omega_reduce();
  if (x.size() != y.size())
    return false;

  // FIXME: could we speed up also using Disj::hash?
  std::list<const Disj*> ptrs;
  for (const auto& dy : y)
    ptrs.push_back(&dy);

  for (const auto& dx : x) {
    auto res = std::find_if(ptrs.begin(), ptrs.end(),
                            [&dx](const Disj* dy) { return dx.equals(*dy); });
    if (res == ptrs.end())
      return false;
    ptrs.erase(res);
  }
  return true;
}

template <typename PH>
bool
PolySet<PH>::geom_covers(const PolySet& y) const {
  const auto& x_seq = seq();
  return all_of(y.seq(),
                [&x_seq](const Disj& d) {
                  return check_containment(d, x_seq);
                });
}

template <typename PH>
dim_type
PolySet<PH>::affine_dim() const {
  if (has_single_disjunct())
    return begin()->affine_dim();

  Disj hull {space_dim(), Spec_Elem::EMPTY, Topol::CLOSED};
  Disj univ {space_dim(), Spec_Elem::UNIVERSE, Topol::CLOSED};

  for (const auto& d : seq()) {
    d.minimize();
    if (d.is_empty())
      continue;
    Disj affine = univ;
    for (const auto& c : d.cons()) {
      if (c.is_equality())
        affine.add_con(c);
    }
    hull.poly_hull_assign(affine);
  }

  return hull.affine_dim();
}

template <typename PH>
Poly_Con_Rel
PolySet<PH>::relation_with(const Con& c) const {
  if (is_empty())
    return Poly_Con_Rel::is_included()
      && Poly_Con_Rel::is_disjoint()
      && Poly_Con_Rel::saturates();

  bool is_included = true;
  bool is_disjoint = true;
  bool is_strictly_intersecting = false;
  bool included_once = false;
  bool disjoint_once = false;
  bool saturates = true;

  for (const auto& d : seq()) {
    const auto rel_i = d.relation_with(c);
    if (rel_i.implies(Poly_Con_Rel::is_included()))
      included_once = true;
    else
      is_included = false;
    if (rel_i.implies(Poly_Con_Rel::is_disjoint()))
      disjoint_once = true;
    else
      is_disjoint = false;
    if (rel_i.implies(Poly_Con_Rel::strictly_intersects()))
      is_strictly_intersecting = true;
    if (!rel_i.implies(Poly_Con_Rel::saturates()))
      saturates = false;
  }

  auto res = Poly_Con_Rel::nothing();
  if (is_included)
    res = res && Poly_Con_Rel::is_included();
  if (is_disjoint)
    res = res && Poly_Con_Rel::is_disjoint();
  if (is_strictly_intersecting || (included_once && disjoint_once))
    res = res && Poly_Con_Rel::strictly_intersects();
  if (saturates)
    res = res && Poly_Con_Rel::saturates();
  return res;
}

template <typename PH>
Poly_Gen_Rel
PolySet<PH>::relation_with(const Gen& g) const {
  for (const auto& d : seq()) {
    const auto rel_i = d.relation_with(g);
    if (rel_i.implies(Poly_Gen_Rel::subsumes()))
      return Poly_Gen_Rel::subsumes();
  }
  return Poly_Gen_Rel::nothing();
}

template <typename PH>
bool
PolySet<PH>::min(const Affine_Expr& ae, Rational& value,
          bool* included_ptr, Gen* g_ptr) const {
  omega_reduce();
  if (is_empty())
    return false;

  // Best value, included and generator.
  Rational best_value;
  bool best_included = false;
  Gen best_g = point();

  // Disjunct-specific value, included and generator.
  Rational d_value;
  bool d_included = false;
  Gen d_g = point();
  bool* d_included_ptr = included_ptr ? &d_included : nullptr;
  Gen* d_g_ptr = g_ptr ? &d_g : nullptr;

  auto update_best = [&]() {
    best_value = std::move(d_value);
    if (included_ptr)
      best_included = d_included;
    if (g_ptr)
      best_g = std::move(d_g);
  };

  bool first = true;
  for (auto& d : seq()) {
    if (not d.min(ae, d_value, d_included_ptr, d_g_ptr))
      return false;
    if (first) {
      first = false;
      update_best();
    } else {
      // not first disjunct: need to compare
      switch (compare(best_value, d_value)) {
      case 1: // found smaller value: update best
        update_best();
        break;
      case 0: // found same value
        if (included_ptr && not best_included && d_included)
          best_included = true;
        if (g_ptr && best_g.is_closure_point() && d_g.is_point())
          best_g = std::move(d_g);
        break;
      case -1: // found bigger value: do nothing
      default:
        break;
      }
    }
  }
  // Move best values found into output parameters.
  value = std::move(best_value);
  if (included_ptr)
    *included_ptr = best_included;
  if (g_ptr)
    *g_ptr = std::move(best_g);
  return true;
}

template <typename PH>
void
PolySet<PH>::set_topology(Topol t) {
  assert(t == Topol::NNC || is_topologically_closed());
  if (impl().topol == t)
    return;
  impl().topol = t;
  std::for_each(seq().begin(), seq().end(),
                [t](Disj& d) { d.set_topology(t); });
  assert(check_inv());
}

template <typename PH>
void
PolySet<PH>::topological_closure_assign() {
  if (is_necessarily_closed())
    return;
  for (auto& d : seq()) {
    if (d.is_topologically_closed())
      continue;
    d.topological_closure_assign();
    clear_reduced();
  }
  assert(check_inv());
}

namespace {

// Note: here the goal is to *efficiently* check if two abstract elements
// are disjoint; we provide two (overloaded) helpers:
// the first one is for *boxed* polyhedra and answers by checking
// the (already cached) bounding box; the second one is meant to
// be used by other abstractions and trivially returns false.

template <typename PH>
inline bool
cheap_is_disjoint_from(const B_Wrap<PH>& x, const B_Wrap<PH>& y) {
  assert(x.has_valid_bbox() && y.has_valid_bbox());
  const auto& x_bb = x.impl_bbox();
  const auto& y_bb = y.impl_bbox();
  return x_bb.is_disjoint_from(y_bb);
}

template <typename PH>
inline bool
cheap_is_disjoint_from(const PH&, const PH&) { return false; }

} // namespace

template <typename PH>
void
PolySet<PH>::intersection_assign(const PolySet& y) {
  auto& x = *this;
  if (y.is_empty()) {
    x.set_empty();
    return;
  }

  x.omega_reduce();
  y.omega_reduce();

  using Ptr = const Disj*;
  using Ptrs = std::vector<Ptr>;

  auto get_seq_pointers = [](const Seq& seq) {
    Ptrs res;
    res.reserve(seq.size());
    for (const auto& d : seq)
      res.push_back(&d);
    return res;
  };

  auto x_ptrs = get_seq_pointers(x.seq());
  auto y_ptrs = get_seq_pointers(y.seq());

  // Note: y_ptrs may contain null pointers.
  auto find_container = [](Ptr x_ptr, Ptrs& y_ptrs) {
    return std::find_if(y_ptrs.begin(), y_ptrs.end(),
                        [x_ptr](auto y_ptr) {
                          return y_ptr != nullptr
                            && y_ptr->contains(*x_ptr);
                        });
  };

  Seq new_seq;

  // Scan x_ptrs and move into new_seq those disjuncts
  // that are covered by a disjunct in y_ptrs.
  for (auto& x_ptr : x_ptrs) {
    auto y_it = find_container(x_ptr, y_ptrs);
    if (y_it == y_ptrs.end())
      continue;
    auto& y_ptr = *y_it;
    // *y_ptr contains *x_ptr : move latter in new_seq
    if (x_ptr->contains(*y_ptr)) {
      // *x_ptr == *y_ptr : avoid rechecking y_ptr
      y_ptr = nullptr;
    }
    new_seq.emplace_back(std::move(*x_ptr));
    // avoid rechecking x_ptr
    x_ptr = nullptr;
  }

  // Scan y_ptrs and move into new_seq those disjuncts
  // that are covered by a disjunct in x_ptrs.
  for (auto& y_ptr : y_ptrs) {
    if (y_ptr == nullptr)
      continue;
    auto x_it = find_container(y_ptr, x_ptrs);
    if (x_it == x_ptrs.end())
      continue;
    assert(*x_it != nullptr && not y_ptr->contains(**x_it));
    // *x_ptr contains *y_ptr : move latter in new_seq
    new_seq.emplace_back(std::move(*y_ptr));
    // avoid rechecking y_ptr
    y_ptr = nullptr;
  }

  // Process remaining elements pairwise.
  for (auto x_ptr : x_ptrs) {
    if (x_ptr == nullptr)
      continue;
    for (auto y_ptr : y_ptrs) {
      if (y_ptr == nullptr)
        continue;
      if (cheap_is_disjoint_from(*x_ptr, *y_ptr))
        continue;
      // Note: copy is meant.
      auto x_d = *x_ptr;
      x_d.intersection_assign(*y_ptr);
      if (x_d.is_empty())
        continue;
      new_seq.emplace_back(std::move(x_d));
    }
  }

  x.seq() = std::move(new_seq);
  x.clear_reduced();
  assert(check_inv());
}

template <typename PH>
PolySet<PH>
PolySet<PH>::split_aux_ineq(const Con& c, Topol t, bool integral) {
  PolySet res(space_dim(), topology(), Spec_Elem::EMPTY);
  if (c.is_tautological())
    return res;
  if (c.is_inconsistent()) {
    m_swap(res);
    return res;
  }
  auto it = seq().begin();
  while (it != seq().end()) {
    auto res_it = integral ? it->integral_split(c) : it->split(c, t);
    res.add_disjunct(std::move(res_it));
    if (it->is_empty())
      it = seq().erase(it);
    else
      ++it;
  }
  assert(check_inv() && res.check_inv());
  return res;
}

template <typename PH>
PolySet<PH>
PolySet<PH>::integral_split_aux_eq(const Con& c) {
  assert(is_necessarily_closed());
  assert(c.is_equality());
  PolySet res(space_dim(), topology(), Spec_Elem::EMPTY);
  if (c.is_tautological())
    return res;
  if (detail::is_integral_inconsistent(c)) {
    m_swap(res);
    return res;
  }
  auto it = seq().begin();
  while (it != seq().end()) {
    // FIXME: ad-hoc 3-way split?
    auto [c_lt, c_gt] = detail::integral_complement_eq(c);
    auto& disj = *it;
    auto res_ge = disj.integral_split(c_lt);
    auto res_lt = std::move(disj);
    auto res_eq = res_ge.integral_split(c_gt);
    auto res_gt = std::move(res_ge);
    res.add_disjunct(std::move(res_lt));
    res.add_disjunct(std::move(res_gt));
    if (res_eq.is_empty())
      it = seq().erase(it);
    else {
      *it = std::move(res_eq);
      ++it;
    }
  }
  assert(check_inv() && res.check_inv());
  return res;
}

template <typename PH>
void
PolySet<PH>::join_assign(const PolySet& y) {
  omega_reduce();
  y.omega_reduce();
  for (const auto& d : y) {
    bool d_known_redundant = false;
    bool d_known_not_redundant = false;
    for (auto it = begin(); it != end(); ) {
      if (not d_known_not_redundant && it->contains(d)) {
        d_known_redundant = true;
        break;
      }
      if (d.contains(*it)) {
        it = seq().erase(it);
        d_known_not_redundant = true;
      } else
        ++it;
    }
    if (not d_known_redundant)
      seq().emplace_back(d);
  }
  assert(check_inv());
}

template <typename PH>
void
PolySet<PH>::remove_higher_space_dims(dim_type new_dim) {
  if (new_dim >= space_dim())
    return;
  for (auto& d : seq()) {
    d.remove_higher_space_dims(new_dim);
    clear_reduced();
  }
  impl().dim = new_dim;
  assert(check_inv());
}

template <typename PH>
void
PolySet<PH>::map_space_dims(const Dims& pfunc) {
  assert(space_dim() == num_rows(pfunc));
  if (space_dim() == 0)
    return;
  auto max_dim = *std::max_element(pfunc.begin(), pfunc.end());
  if (max_dim == not_a_dim()) {
    impl().dim = 0;
    if (is_empty())
      set_empty();
    else
      set_universe();
    return;
  }

  const auto old_dim = space_dim();
  impl().dim = max_dim + 1;
  for (auto& disj : seq())
    disj.map_space_dims(pfunc);
  if (space_dim() < old_dim)
    clear_reduced();
  assert(check_inv());
}

template <typename PH>
void
PolySet<PH>::widening_assign(const PolySet& y, const Cons* upto_ptr,
                             Widen_Impl w_impl, Widen_Spec w_spec) {
  auto& x = *this;
  if (y.is_empty())
    return;

  if (w_spec == Widen_Spec::SAFE) {
    // Trivial lifting of risky widening.
    x.join_assign(y);
  }

  auto x_hull = Disj{x.space_dim(), x.topology(), Spec_Elem::EMPTY};
  for (const auto& x_disj : x)
    x_hull.poly_hull_assign(x_disj);

  auto y_hull = Disj{y.space_dim(), y.topology(), Spec_Elem::EMPTY};
  for (const auto& y_disj : y)
    y_hull.poly_hull_assign(y_disj);

  // first case
  using Cert = bhrz03_widen::Cert;
  const Cert x_cert {x_hull};
  const Cert y_cert {y_hull};
  int convergence = y_cert.compare(x_cert);

  // check if certificate is satisfied
  if (convergence == 1 ||
     (convergence == 0 && y.size() > 1 && x.size() == 1))
    return;
  else if ((convergence == 0) &&
           (x.size() > 1 && y.size() > 1)) {
    auto cert_comp
      = [] (const Cert& x_cert, const Cert& y_cert) -> bool {
          return (x_cert.compare(y_cert) == 1);
        };

    using Cert_ms = std::map<Cert, std::size_t, decltype(cert_comp)>;
    auto collect_certs
      = [](const PolySet& ps, Cert_ms& cert_ms) {
          for (const auto& d : ps) {
            // if (d.impl().cs.sk_rows.size() <= 1)
            //   continue;
            auto ret_val = cert_ms.emplace(Cert{d}, std::size_t{});
            ++(ret_val.first->second);
          }
        };

    Cert_ms x_cert_ms(cert_comp);
    Cert_ms y_cert_ms(cert_comp);

    collect_certs(x, x_cert_ms);
    collect_certs(y, y_cert_ms);

    auto test_convergence
      = [](const Cert_ms& x_ms, const Cert_ms& y_ms) -> bool {
          auto xi = x_ms.begin(), x_last = x_ms.end(),
            yi = y_ms.begin(), y_last = y_ms.end();
          while (xi != x_last && yi != y_last) {
            const Cert& x_cert = xi->first;
            const Cert& y_cert = yi->first;
            switch (x_cert.compare(y_cert)) {
            case 0:
              {
                const size_t& x_count = xi->second;
                const size_t& y_count = yi->second;
                if (x_count == y_count) {
                  ++xi;
                  ++yi;
                }
                else
                  return (x_count < y_count);
                break;
              }
            case 1:
              return false;
            case -1:
              return true;
            }
          }
          return yi != y_last;
        };

    if (test_convergence(x_cert_ms, y_cert_ms))
      return;
  }

  // second case
  if (x_hull.strictly_contains(y_hull)) {
    auto x_hull_old = x_hull;
    if (upto_ptr == nullptr)
      x_hull.widening_assign(y_hull, w_impl, Widen_Spec::RISKY);
    else
      x_hull.widening_assign(y_hull, *upto_ptr, w_impl, Widen_Spec::RISKY);
    x_hull.poly_difference_assign(x_hull_old);
    PolySet d {x_hull};
    x.join_assign(d);
    return;
  }

  // third case
  PolySet new_x {x_hull};
  x = new_x;
  return;
}

template <typename PH>
void
PolySet<PH>::collapse(dim_type n) {
  // If n <= 0, do nothing at all.
  if (n <= 0)
    return;
  // Try to be lazy with (empty/omega) reductions.
  if (num_disjuncts() <= n)
    return;
  empty_reduce();
  if (num_disjuncts() <= n)
    return;
  omega_reduce();
  if (num_disjuncts() <= n)
    return;

  // FIXME: find a reasonable criterion to guide merging of disjuncts.
  // For the time being, we simply collapse the trailing part of
  // the list of disjuncts.
  auto nth = seq().begin();
  std::advance(nth, n-1);
  auto last = seq().end();
  assert(nth != last);
  auto next = nth;
  ++next;
  if (next == last)
    return;
  for (auto it = next; it != last; ++it)
    nth->poly_hull_assign(*it);
  seq().erase(next, last);
  if (n == 1)
    set_reduced();
  else
    clear_reduced();
  assert(check_inv());
}

template <typename PH>
void
PolySet<PH>::add_space_dims(dim_type m, bool project) {
  assert(m >= 0);
  if (m == 0)
    return;
  if (is_empty()) {
    impl().dim += m;
    return;
  }
  for (auto& d : seq())
    d.add_space_dims(m, project);
  impl().dim += m;
  assert(check_inv());
}

template <typename PH>
void
PolySet<PH>::difference_assign(const PolySet& y) {
  auto& x = *this;
  using std::swap;
  x.omega_reduce();
  y.omega_reduce();
  Seq new_seq = x.seq();
  for (const auto& yd : y.seq()) {
    Seq tmp_seq;
    for (const auto& nd : new_seq) {
      auto p = linear_partition(yd, nd);
      tmp_seq.splice(tmp_seq.end(), p.second);
    }
    swap(tmp_seq, new_seq);
  }
  swap(x.seq(), new_seq);
  x.clear_reduced();
  assert(x.check_inv());
}

template <typename PH>
void
PolySet<PH>::concatenate_assign(const PolySet& y) {
  if (is_empty()) {
    impl().dim += y.space_dim();
    return;
  }
  if (y.is_empty()) {
    impl().dim += y.space_dim();
    set_empty();
    return;
  }
  const auto& x = *this;
  x.omega_reduce();
  y.omega_reduce();
  PolySet z(x.space_dim() + y.space_dim(), x.topology(), Spec_Elem::EMPTY);
  for (const auto& xd : x) {
    for (const auto& yd : y) {
      z.seq().emplace_back(xd);
      z.seq().back().concatenate_assign(yd);
    }
  }
  z.set_reduced();
  m_swap(z);
  assert(check_inv());
}

template <typename PH>
void
PolySet<PH>::expand_space_dim(Var var, dim_type m) {
  for (auto& disj : seq())
    disj.expand_space_dim(var, m);
  impl().dim += m;
  assert(check_inv());
}

template <typename PH>
void
PolySet<PH>::fold_space_dims(const Index_Set& vars, Var dest) {
  auto num_folded = vars.size();
  if (num_folded == 0)
    return;
  for (auto& d : seq())
    d.fold_space_dims(vars, dest);
  impl().dim -= num_folded;
  clear_reduced();
  assert(check_inv());
}

template <typename PH>
bool
PolySet<PH>::ascii_load(std::istream& s) {
  auto sz = 0;
  if (not (ascii_load_string(s, "size") && (s >> sz)))
    return false;

  auto sd = 0;
  if (not (ascii_load_string(s, "dim") && (s >> sd)))
    return false;

  auto topol = Topol::CLOSED;
  std::string str;
  if (not (ascii_load_string(s, "topol") && (s >> str)))
    return false;
  if (str == "CLOSED")
    topol = Topol::CLOSED;
  else if (str == "NNC")
    topol = Topol::NNC;
  else
    return false;

  auto red = 0;
  if (not (ascii_load_string(s, "reduced") && (s >> red)))
    return false;
  if (not (red == 0 || red == 1))
    return false;

  PolySet new_set {sd, topol, Spec_Elem::EMPTY};
  Disj d;
  for (auto i = 0; i != sz; ++i) {
    auto idx = 0;
    if (not (ascii_load_string(s, "===") &&
             ascii_load_string(s, "start") &&
             ascii_load_string(s, "disjunct") &&
             (s >> idx) && (idx == i) &&
             ascii_load_string(s, "===")))
      return false;
    if (!d.ascii_load(s))
      return false;
    new_set.seq().push_back(std::move(d));
    if (not (ascii_load_string(s, "===") &&
             ascii_load_string(s, "end") &&
             ascii_load_string(s, "disjunct") &&
             (s >> idx) && (idx == i) &&
             ascii_load_string(s, "===")))
      return false;
  }
  new_set.impl_.reduced = static_cast<bool>(red);
  swap(*this, new_set);
  assert(check_inv());
  return true;
}

template <typename PH>
void
PolySet<PH>::ascii_dump(std::ostream& s) const {
  using namespace IO_Operators;
  s << "size " << size() << "\n";
  s << "dim " << space_dim() << "\n";
  s << "topol " << topology() << "\n";
  s << "reduced " << is_reduced() << "\n\n";

  auto i = 0;
  for (const auto& d : seq()) {
    s << "=== start disjunct << " << i << " ===\n";
    d.ascii_dump(s);
    s << "=== end disjunct << " << i << " ===\n\n";
    ++i;
  }
}

// Explicit template instantiation (definition).
template class PolySet<B_Poly>;
template class PolySet<BF_Poly>;

} // namespace pplite
