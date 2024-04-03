/* PPLite: a lightweight library for convex polyhedra derived from PPL.
   Copyright (C) 2001-2010 Roberto Bagnara <bagnara@cs.unipr.it>
   Copyright (C) 2010-2018 BUGSENG srl (http://bugseng.com)
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

#ifndef pplite_Poly_templ_hh
#define pplite_Poly_templ_hh

#include "Poly.hh"
#include "Poly_min.hh"

#include <algorithm>
#include <functional>
#include <numeric> // for std::iota
#include <type_traits>

namespace pplite {

namespace detail {

template <typename Dst, typename Src>
void
concat_sys(Dst& dst, Src&& src) {
  // Add sing elements.
  std::move(src.sing_rows.begin(), src.sing_rows.end(),
            std::back_inserter(dst.sing_rows));
  // Add sk elements.
  const dim_type old_sk_size = dst.sk_rows.size();
  std::move(src.sk_rows.begin(), src.sk_rows.end(),
            std::back_inserter(dst.sk_rows));
  // Add ns elements, shifting their indices.
  auto pos = num_rows(dst.ns_rows);
  for (auto i : index_range(src.ns_rows)){
    dst.ns_rows.push_back(std::move(src.ns_rows[i]));
    dst.ns_rows[pos++] >>= old_sk_size;
  }
}

template <typename Row>
inline void
permute_space_dims_cycle(std::vector<Row>& rows,
                         const Dims& cycle, dim_type d) {
  assert(cycle.size() > 1);
  assert(d >= 1 + *(std::max_element(cycle.begin(), cycle.end())));
  for (auto& r : rows)
    r.permute_space_dims_cycle(cycle, d);
}

template <typename Rows>
inline void
permute_space_dims_cycle(Poly_Impl::Sys<Rows>& sys,
                         const Dims& cycle, dim_type d) {
  permute_space_dims_cycle(sys.sing_rows, cycle, d);
  permute_space_dims_cycle(sys.sk_rows, cycle, d);
  // ns_rows are unaffected.
}

inline void
permute_space_dims_cycle(Poly_Impl& ph, const Dims& cycle, dim_type d) {
  permute_space_dims_cycle(ph.cs, cycle, d);
  permute_space_dims_cycle(ph.cs_pending, cycle, d);
  permute_space_dims_cycle(ph.gs, cycle, d);
  permute_space_dims_cycle(ph.gs_pending, cycle, d);
  // sat matrices are unaffected.
}

template <typename T>
void
permute_space_dims(T& t, const Dims& perm, dim_type d) {
  assert(d >= num_rows(perm));
  Dims cycle;
  cycle.reserve(d);
  std::vector<bool> visited(d);
  for (auto i : bwd_range(d)) {
    if (visited[i])
      continue;
    dim_type j = i;
    do {
      visited[j] = true;
      const dim_type k = perm[j];
      assert(k != not_a_dim());
      if (k == j)
        break;
      cycle.push_back(j);
      j = k;
    } while (!visited[j]);
    if (cycle.size() > 1)
      permute_space_dims_cycle(t, cycle, d);
    cycle.clear();
  }
}

inline void
remove_invalid_lines(Gens& sing_rows) {
  Index_Set tbr = invalid_lines(sing_rows);
  erase_using_sorted_indices(sing_rows, tbr);
}

inline void
remove_invalid_rays(Gens& sk_rows, NS_Rows& ns_rows) {
  Index_Set tbr = invalid_rays(sk_rows);
  if (tbr.empty())
    return;
  erase_using_sorted_indices(sk_rows, tbr);
  for (auto& ns : ns_rows)
    ns.remove_all(tbr);
  promote_singletons(ns_rows, sk_rows);
}

inline void
remove_invalid_lines_and_rays(Poly_Impl::Sys<Gens>& gs) {
  remove_invalid_lines(gs.sing_rows);
  remove_invalid_rays(gs.sk_rows, gs.ns_rows);
}

template <typename Rows>
inline void
erase_higher_dims(Poly_Impl::Sys<Rows>& sys, dim_type new_dim) {
  erase_higher_dims(sys.sing_rows, new_dim);
  erase_higher_dims(sys.sk_rows, new_dim);
}

template <typename Rows, typename Iter>
inline void
erase_space_dims(Poly_Impl::Sys<Rows>& sys, Iter first, Iter last) {
  erase_space_dims(sys.sing_rows, first, last);
  erase_space_dims(sys.sk_rows, first, last);
}

template <typename PH>
void
par_affine_image_aux(PH& ph,
                     const Vars& vars, const Linear_Exprs& exprs,
                     const Integers& inhomos, const Integers& dens) {
  // This implementation starts by computing a minimized DD pair,
  // since it is meant to avoid non-incremental calls to the
  // conversion algorithm.
  ph.minimize();
  if (ph.is_empty())
    return;

  const auto n_vars = num_rows(vars);
  if (n_vars == 1) {
    ph.affine_image(vars[0], exprs[0], inhomos[0], dens[0]);
    return;
  }

  // Terminology: for each b_i in [0, n_vars), we say that
  //   (vars[b_i], exprs[b_i], inhomos[b_i], dens[b_i])
  // is a *binding*, having *lhs* vars[b_i].id() and *rhs* exprs[b_i].

  // Build the "binding index" dependency graph: there is an arc
  // b_1 --> b_2 if the rhs of b_2 mentions the lhs of b_1;
  // intuitively, b_2 depends on the old value of the lhs of b_1,
  // hence it has to be evaluated before evaluating b_1.
  // For each index b_i, deps[b_i] stores the b_j's such that b_i --> b_j.
  // If deps[b_i] is empty, the binding is free from dependencies
  // and thus it can be evaluated.

  // The graph of dependencies.
  using Deps = std::vector<Index_Set>;
  Deps deps;
  // The bindings still waiting to be processed.
  Index_Set waiting;

  // Direct and inverse maps from binding index to lhs.
  const auto sdim = ph.space_dim();
  Dims lhs2bnd(sdim, not_a_dim());
  for (auto i : bwd_range(n_vars))
    lhs2bnd[vars[i].id()] = i;
  auto binding_to_lhs = [&vars](dim_type i) { return vars[i].id(); };
  auto lhs_to_binding = [&lhs2bnd](dim_type i) { return lhs2bnd[i]; };

  // Helper lambda to build the graph.
  auto build_graph = [&]() mutable {
    deps = Deps(n_vars);
    for (auto b_i : range(n_vars)) {
      const auto& rhs_i = exprs[b_i];
      for (auto lhs_j : bwd_dim_range(rhs_i)) {
        auto b_j = lhs_to_binding(lhs_j);
        // Skip those rhs vars that do not occur in lhs.
        if (b_j == not_a_dim())
          continue;
        if (rhs_i.get(lhs_j) == 0)
          continue;
        // b_i depends on b_j
        deps[b_j].set(b_i);
      }
    }
  };

  // Helper lambdas to classify bindings.
  auto is_invertible = [&](dim_type b) { return deps[b].test(b); };
  auto is_identity = [&](dim_type b) {
    auto lhs = binding_to_lhs(b);
    const auto& rhs = exprs[b];
    return is_invertible(b)
      && rhs.get(lhs) == 1 && inhomos[b] == 0 && dens[b] == 1
      && rhs.all_zeroes(0, lhs)
      && rhs.all_zeroes(lhs + 1, rhs.space_dim());
  };
  auto is_processable_noninv = [&](dim_type b) {
    return deps[b].empty();
  };
  auto is_processable_inv = [&](dim_type b) {
    return is_invertible(b) && deps[b].size() == 1;
  };

  // Helper lambdas to update graph after processing a binding.
  auto remove_from_deps = [&](dim_type b) mutable {
    // Note: we only update the deps for unprocessed bindings.
    for (auto bj : waiting)
      deps[bj].reset(b);
  };
  auto remove_all_from_deps = [&](const Index_Set& bs) mutable {
    // Note: we only update the deps for unprocessed bindings.
    for (auto bj : waiting)
      deps[bj] -= bs;
  };

  // Helper lambda to process bindings.
  auto process_noninv = [&](const Index_Set& bs) mutable {
    assert(all_of(bs, is_processable_noninv));
    for (auto b : bs)
      ph.unconstrain(vars[b]);
    for (auto b : bs)
      ph.add_con(dens[b] * vars[b] - exprs[b] == inhomos[b]);
    // Update waiting set and dependency graph.
    waiting -= bs;
    remove_all_from_deps(bs);
  };
  auto process_inv = [&](const Index_Set& bs) mutable {
    assert(all_of(bs, is_processable_inv));
    for (auto b : bs) {
      ph.affine_image(vars[b], exprs[b], inhomos[b], dens[b]);
      // Update waiting set and dependency graph.
      waiting.reset(b);
      remove_from_deps(b);
    }
  };

  // Init the dependency graph.
  build_graph();
  // Init the waiting set (ignore identities).
  for (auto b_i : bwd_range(n_vars)) {
    if (!is_identity(b_i))
      waiting.set(b_i);
  }
  // Process the bindings still waiting.
  Dims added_dims;
  while (!waiting.empty()) {
    // Collect all the processable bindings.
    Index_Set proc_noninv;
    Index_Set proc_inv;
    for (auto w_b : waiting) {
      if (is_processable_noninv(w_b))
        proc_noninv.set(w_b);
      else if (is_processable_inv(w_b))
        proc_inv.set(w_b);
    }
    // Process them.
    if (!(proc_noninv.empty() && proc_inv.empty())) {
      process_noninv(proc_noninv);
      process_inv(proc_inv);
      continue;
    }
    // No binding is processable (i.e., graph has a non-trivial cycle).
    assert(!waiting.empty());
    // Break a cycle by processing a binding into a new space dimension
    // (prefer the binding having most dependencies).
    auto b = *std::max_element(waiting.begin(), waiting.end(),
                               [&](dim_type bi, dim_type bj) {
                                 return deps[bi].size() < deps[bj].size();
                               });
    // Add a new space dim (for the lhs of b)
    // and process the rewritten version of binding b.
    ph.add_space_dims(1);
    Var new_var(sdim + num_rows(added_dims));
    ph.add_con(dens[b] * new_var - exprs[b] == inhomos[b]);
    // Remember the identity of the lhs of b (causing new space dim).
    added_dims.push_back(vars[b].id());
    // Update waiting and dependency graph.
    waiting.reset(b);
    remove_from_deps(b);
  }
  // Now remap added dims, if any.
  if (added_dims.empty())
    return;
  auto new_sdim = ph.space_dim();
  assert(new_sdim == sdim + num_rows(added_dims));
  Dims pfunc(new_sdim);
  std::iota(pfunc.begin(), pfunc.end(), 0);
  for (auto i : index_range(added_dims)) {
    pfunc[sdim + i] = added_dims[i];
    pfunc[added_dims[i]] = not_a_dim();
  }
  ph.map_space_dims(pfunc);
  assert(ph.space_dim() == sdim);
}

} // namespace detail

template <typename Iter>
void
Poly_Impl::remove_space_dims(Iter first, Iter last) {
  using value_type = typename std::iterator_traits<Iter>::value_type;
  static_assert(std::is_same<value_type, dim_type>::value,
                "Iter must have dim_type as its value type");
  assert(std::all_of(first, last,
                     [this](value_type i) { return 0 <= i && i < dim; }));
  assert(std::is_sorted(first, last));
  const dim_type num_rem = std::distance(first, last);
  if (num_rem == 0)
    return;
  const dim_type new_dim = dim - num_rem;
  ensure_valid_gens();
  if (marked_empty()) {
    dim = new_dim;
    return;
  }
  if (new_dim == 0) {
    dim = 0;
    set_universe();
    return;
  }

  const bool removing_higher = (*first == new_dim);

  const auto rem_percentage = (num_rem * 100 / dim);
  const bool preserve_dd = (rem_percentage < remove_space_dims_percentage);

  if (preserve_dd) {
    unconstrain(first, last);
    minimize();
    if (removing_higher) {
      detail::erase_higher_dims(cs, new_dim);
      detail::erase_higher_dims(gs, new_dim);
    } else {
      detail::erase_space_dims(cs, first, last);
      detail::erase_space_dims(gs, first, last);
    }
    // Remove invalid lines (no need to check rays).
    detail::remove_invalid_lines(gs.sing_rows);
    dim = new_dim;
  } else {
    // Only preserve the generators.
    Sys<Gens> tmp_gs = std::move(gs);
    detail::concat_sys(tmp_gs, std::move(gs_pending));
    if (removing_higher)
      detail::erase_higher_dims(tmp_gs, new_dim);
    else
      detail::erase_space_dims(tmp_gs, first, last);
    // Remove invalid generators.
    detail::remove_invalid_lines_and_rays(tmp_gs);
    // Use tmp to reinit the polyhedron.
    reinit_with_gens(new_dim, topol, std::move(tmp_gs));
  }

  assert(check_inv());
}

template <typename Iter>
void
Poly_Impl::add_cons(Iter first, Iter last) {
  assert(std::all_of(first, last,
                     [this](const Con& c) {
                       return c.space_dim() <= this->space_dim();
                     }));
  if (marked_empty())
    return;
  if (space_dim() == 0) {
    if (std::any_of(first, last, std::mem_fn(&Con::is_inconsistent)))
      set_empty();
    return;
  }

#ifndef NDEBUG
  auto is_non_trivial_strict = [](const Con& c) {
    return c.is_strict_inequality()
           && !c.is_tautological() && !c.is_inconsistent();
  };
  assert(!is_necessarily_closed() ||
         std::none_of(first, last, is_non_trivial_strict));
#endif

  ensure_valid_cons();
  auto& eq = cs_pending.sing_rows;
  auto& sk = cs_pending.sk_rows;
  for ( ; first != last; ++first) {
    const auto& c = *first;
    if (c.is_equality())
      eq.push_back(c);
    else
      sk.push_back(c);
  }
  if (marked_min() && !cs_pending.empty())
    set_status(Status::PENDING);
  assert(check_inv());
}

template <typename Iter>
void
Poly_Impl::add_gens(Iter first, Iter last) {
  assert(std::all_of(first, last,
                     [this](const Gen& g) {
                       return g.space_dim() <= this->space_dim();
                     }));
  assert(!is_necessarily_closed() ||
         std::none_of(first, last, std::mem_fn(&Gen::is_closure_point)));

  ensure_valid_gens();
  if (marked_empty()) {
    auto iter = std::find_if(first, last, std::mem_fn(&Gen::is_point));
    assert(iter != last);
    add_gen(*iter);
    add_gens(first, iter);
    ++iter;
    add_gens(iter, last);
    return;
  }
  if (dim == 0)
    return;
  auto& lines = gs_pending.sing_rows;
  auto& sk = gs_pending.sk_rows;
  for ( ; first != last; ++first) {
    const auto& g = *first;
    if (g.is_line())
      lines.push_back(g);
    else
      sk.push_back(g);
  }
  if (marked_min() && !gs_pending.empty())
    set_status(Status::PENDING);
  assert(check_inv());
}

template <typename System, typename Out, typename Pred>
Out
Poly_Impl::materialize(const System& sys, const System& pending, Out out,
                       bool cons_and_closed, Pred pred) {
  // Copy singular.
  out = std::copy(sys.sing_rows.begin(), sys.sing_rows.end(), out);
  out = std::copy(pending.sing_rows.begin(), pending.sing_rows.end(), out);
  if (cons_and_closed) {
    // Materializing constraints for a CLOSED poly.
    // Copy skel, but not the (strict) positivity.
    out = std::copy_if(sys.sk_rows.begin(), sys.sk_rows.end(),
                       out, pred);
    out = std::copy_if(pending.sk_rows.begin(), pending.sk_rows.end(),
                       out, pred);
    // No need to materialize non-skel.
    return out;
  }
  // Materializing generators of constraints for an NNC poly.
  // Copy skel.
  out = std::copy(sys.sk_rows.begin(), sys.sk_rows.end(), out);
  out = std::copy(pending.sk_rows.begin(), pending.sk_rows.end(), out);
  // Materialize non-skel.
  for (const auto& ns : sys.ns_rows) {
    *out = detail::materialize(ns, sys.sk_rows);
    ++out;
  }
  for (const auto& ns : pending.ns_rows) {
    *out = detail::materialize(ns, pending.sk_rows);
    ++out;
  }
  return out;
}

template <typename Iter>
void
Poly_Impl::unconstrain(Iter first, Iter last) {
  assert(std::all_of(first, last,
                     [this](dim_type d) { return d < space_dim(); }));
  if (first == last)
    return;
  ensure_valid_gens();
  if (marked_empty())
    return;
  auto& lines = gs_pending.sing_rows;
  for (auto it = first; it != last; ++it)
    lines.push_back(line(Var(*it)));
  if (marked_min())
    set_status(Status::PENDING);
  assert(check_inv());
}

/**********************************************************************/

namespace detail {

/* Structures needed for n-ary constraint hull */

// Constraint with a normalized slope.
struct Norm_Con {
  Linear_Expr slope;
  Rational inhomo;

  // Indices of the polys against which it has already been tested.
  Index_Set valid_polys;

  Norm_Con(const Con& c) : slope(c.linear_expr()) {
    Integer dummy = 0;
    slope.normalize(dummy);
    const Integer& num = c.inhomo_term();
    Integer den = c.linear_expr().gcd(0, slope.space_dim());
    inhomo = Rational(num, den);
  }

  explicit Norm_Con(Linear_Expr&& s)
    : slope(std::move(s)), inhomo(Rational(0)) { }

  Con get_con() const {
    return (inhomo.get_den() * slope + inhomo.get_num() >= 0);
  }
};

using Norm_Cons = std::vector<Norm_Con>;

// Computes the cons of the most precise constraint-hull of a sequence of polys.
template <typename Iter>
Cons
con_hull_cons(Iter first, Iter last, bool boxed) {
  assert(first != last);

  // Collects all the constraints of the polys in [first, last)
  // merging the ones sharing the same slope into their loosest version.
  // For each constraint stores the polys for which it is syntactically valid.
  auto collect_all_cons = [first, last, boxed] () {

    // Helper: merge a constraint with the others that share the same slope.
    auto add_inequality = [] (const Norm_Con& norm_c,
                              Norm_Cons& cons) {
      for (auto& c1 : cons) {
        if (c1.slope.is_equal_to(norm_c.slope)) {

          // Check if it still needs to be initialized.
          if (c1.valid_polys.empty()) {
            c1 = norm_c;
            return;
          }

          // Maybe-update the stored inhomo, loosening the con.
          c1.inhomo = std::max(c1.inhomo, norm_c.inhomo);
          c1.valid_polys |= norm_c.valid_polys;
          return;
        }
      } // end loop on cons

      // Not found. Insert a new slope.
      cons.push_back(norm_c);
    };

    // Helper: add default interval slopes.
    auto add_box_constraints = [first] (Norm_Cons& cs) {
      for (auto d : dim_range(*first)) {
        cs.push_back(Norm_Con(Linear_Expr(Var(d))));
        cs.push_back(Norm_Con(Linear_Expr(-Var(d))));
      }
    };

    /*****************************/

    Norm_Cons cs;

    if (boxed)
      add_box_constraints(cs);

    auto index = 0;
    for (auto it = first; it != last; ++it, ++index) {
      for (const auto& c : it->normalized_cons()) {
        if (c.space_dim() == 0)
          continue;

        Norm_Con norm_c(c);
        norm_c.valid_polys.set(index);

        // Always split equalities.
        // TODO: consider to keep normalized equalities too.

        add_inequality(norm_c, cs);

        if (c.is_equality()) {
          norm_c.slope = - (norm_c.slope);
          norm_c.inhomo = - (norm_c.inhomo);
          add_inequality(norm_c, cs);
        }

      } // end loop on it->cons
    } // end loop on polys

    return cs;
  }; // end of collect cons helper.

  // Updates `c.valid_polys' checking a semantic validity.
  auto store_valid_polys = [first, last] (Norm_Con& c) {
    Con con = c.get_con();
    auto index = 0;
    for (auto it = first; it != last; ++it, ++index) {
      if (c.valid_polys.test(index))
        continue;

      // Waits to be initialized.
      if (c.valid_polys.empty())
        continue;

      if (it->relation_with(con).implies(Poly_Con_Rel::is_included()))
        c.valid_polys.set(index);
    }
  }; // end of store_valid_polys helper.

  // Updates `c' as the stronger non-strict inequality constraint s.t.:
  //  *) has the same slope of c;
  //  *) is not stronger than c;
  //  *) is valid for all polys in [first, last).
  auto make_bounding_con = [first, last] (Norm_Con& c) {
    auto index = 0;
    for (auto it = first; it != last; ++it, ++index) {
      if (c.valid_polys.test(index))
        continue;

      Rational value;
      auto ae = Affine_Expr(c.inhomo.get_den() * c.slope,
                            c.inhomo.get_num());
      bool nonempty_and_bounded = it->min(ae, value);
      if (!nonempty_and_bounded)
        // Cannot be moved.
        return false;

      // Still needs to be initialized.
      if (c.valid_polys.empty()) {
        c.inhomo = -value;
        c.valid_polys.set(index);
        continue;
      }

      // Make sure that we are not strengthening c.
      if (sgn(value) < 0) {
        value /= Rational(c.inhomo.get_den());
        c.inhomo -= value;
      }

      c.valid_polys.set(index);
    }

    // Here it has been successfully updated against all polys.
    return true;
  }; // end of make_bounding_con helper.

  /**********************************************/

  for (auto it = first; it != last; ++it)
    it->minimize();

  // Collect all the constraints checking syntactic validity.
  auto all_cons = collect_all_cons();

  // Check and move the constraints.
  Cons res;
  for (auto& c : all_cons) {
    store_valid_polys(c);
    // Move against invalid polys.
    bool exists = make_bounding_con(c);
    if (exists)
      res.push_back(c.get_con());
  }
  return res;
}

} // namespace detail

template <typename Iter, typename PH>
inline void
con_hull(PH& ph, Iter first, Iter last, bool boxed) {
  if (first == last)
    return;

  auto not_empty = [] (const PH& p) { p.minimize(); return !p.is_empty(); };
  SequenceAdapter<PH> args;
  args.append(first, last, not_empty);
  if (not_empty(ph))
    args.push_back(&ph);

  Cons cs = detail::con_hull_cons(args.begin(), args.end(), boxed);
  ph.set_universe();
  ph.add_cons(std::move(cs));
}

template <typename Iter, typename PH>
inline PH
con_hull(Iter first, Iter last, bool boxed) {
  if (first == last)
    return PH(0, Spec_Elem::EMPTY);

  Cons cs = detail::con_hull_cons(first, last, boxed);
  PH ph(first->space_dim(), Spec_Elem::UNIVERSE, first->topology());
  ph.add_cons(std::move(cs));

  return ph;
}

inline void
Poly::con_hull_assign(const Poly& y, bool boxed) {
  con_hull(*this, &y, &y + 1, boxed);
}

namespace detail {

template <typename PH>
PH
poly_difference(PH& x, const PH& y) {
  if (x.is_empty() || y.is_empty())
    return x;

  PH res(x.space_dim(), Spec_Elem::EMPTY, x.topology());
  if (x.space_dim() == 0 || y.contains(x))
    return res;

  if (x.topology() == Topol::CLOSED) {
    // Here we can't add strict inequalities (when complementing
    // equalities). Thus, if x does not satisfy an equality of y,
    // the result is x.
    for (const auto& c : y.cons()) {
      if (c.is_equality() &&
          !x.relation_with(c).implies(Poly_Con_Rel::is_included()))
        return x;
    }
  }

  for (const auto& c : y.cons()) {
    assert(!c.is_tautological() && !c.is_inconsistent());
    if (x.topology() == Topol::CLOSED && c.is_equality())
      // case already dealt with
      continue;
    if (x.relation_with(c).implies(Poly_Con_Rel::is_included()))
      continue;

    // Build the "complement" of c.
    auto c_neg = detail::complement_con(c, x.topology());

    if (c.is_inequality()) {
      auto z = x;
      z.add_con(c_neg);
      res.poly_hull_assign(z);
    } else {
      // Equality: have to intersect with the two open halfspaces.
      // Speculative optimization: use split + add_con.
      assert(x.topology() == Topol::NNC);
      assert(c_neg.is_strict_inequality());
      auto z1 = x;
      auto z2 = z1.split(c_neg);
      // Build the opposite open halfspace
      auto c_neg2 = Con(-(c_neg.linear_expr()), -(c_neg.inhomo_term()),
                        Con::STRICT_INEQUALITY);
      z2.add_con(c_neg2);
      res.poly_hull_assign(z1);
      res.poly_hull_assign(z2);
    }
  }
  return res;
}

} // namespace detail

inline void
Poly::poly_difference_assign(const Poly& y) {
  auto& x = *this;
  x = detail::poly_difference(x, y);
}

namespace detail {

// Low level helper. Assumptions:
//  - nx in { 0, 1 }
//  - x is not empty and "simple", i.e., defined by nx constraints;
//  - y is not empty.
template <typename PH>
bool
simple_poly_is_disjoint_from_poly(dim_type nx, const PH& x, const PH& y) {
  assert(nx == 0 || nx == 1);
  assert(not x.is_empty() && not y.is_empty());
  if (nx == 0) {
    assert(x.is_universe());
    return false;
  }

  const auto& x_cons = x.cons();
  const auto& x_c = *(x_cons.begin());
  auto x_itv = itv_from_con_inhomo(x_c);
  auto ae = Affine_Expr(x_c.linear_expr());
  auto y_itv = y.get_bounds(ae);
  if (y_itv.is_disjoint_from(x_itv))
    return true;
  // In general, NNC polyhedra may be disjoint even when the
  // itvs are not disjoint (they touch on one boundary).
  // We need to perform a case analysis.
  // 1) polyhedra are *necessarily* closed ==> not disjoint
  if (x.is_necessarily_closed()) {
    assert(y.is_necessarily_closed());
    assert(!x_c.is_strict_inequality());
    return false;
  }
  // 2) both itvs are singletons ==> not disjoint
  if (x_c.is_equality() && y_itv.is_singleton()) {
    assert(x_itv.is_singleton());
    return false;
  }
  // 3) NNC polyhedra are topologically closed ==> not disjoint
  if (not x_c.is_strict_inequality() && y.is_topologically_closed()) {
    assert(x.is_topologically_closed());
    return false;
  }
  // compute itvs intersection in z_itv
  auto z_itv = x_itv;
  z_itv.glb_assign(y_itv);
  assert(not z_itv.is_empty());
  // 4) itvs properly overlap ==> not disjoint
  if (not z_itv.is_singleton())
    return false;
  // 5) c_inhomo_itv and z_itv have different lb ==> not disjoint
  if (x_itv.lb != z_itv.lb)
    return false;
  // 6) x_itv and z_itv have same lb and c is strict ==> disjoint
  if (x_c.is_strict_inequality())
    return true;
  // 7) use min/max to check inclusion in y of the common boundary
  Rational value;
  bool is_bounded = false;
  bool is_included = false;
  if (y_itv.has_lb() && (y_itv.lb == z_itv.lb)) {
    is_bounded = y.min(ae, value, &is_included);
    assert(is_bounded && value == y_itv.lb);
  } else {
    assert(y_itv.has_ub() && (y_itv.ub == z_itv.lb));
    is_bounded = y.max(ae, value, &is_included);
    assert(is_bounded && value == y_itv.ub);
  }
  (void) is_bounded;
  return !is_included;
}

} // namespace detail

} // namespace pplite

#endif // !defined(pplite_Poly_templ_hh)

