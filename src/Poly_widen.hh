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

#ifndef pplite_Poly_widen_hh
#define pplite_Poly_widen_hh 1

#include "globals.hh"
#include "Poly.hh"
#include "Poly_min.hh"
#include "Poly_templ.hh"
#include "B_Poly.hh"
#include "F_Poly.hh"
#include "Scalar_Prod.hh"
#include <map>
#include <vector>

namespace pplite {

namespace detail {

// Check if joining y with x would change (i.e., increase)
// the affine dimension of y; ad-hoc implementation for Poly_Impl.
inline bool
increases_affine_dim(const Poly_Impl& x, const Poly_Impl& y) {
  // Check that all equalities of y are satisfied by all lines
  // and skel gens of x (no need to check non-skel gens).
  const auto& y_eqns = y.cs.sing_rows;
  return not (sp::satisfied_by_all(y_eqns, x.gs.sing_rows)
              && sp::satisfied_by_all(y_eqns, x.gs.sk_rows));
}

// Check if joining y with x would change (i.e., increase)
// the affine dimension of y; generic implementation.
template <typename PH>
inline bool
increases_affine_dim(const PH& x, const PH& y) {
  // Check that all equalities of y are satisfied by all gens of x.
  // TODO: there is no need to scan and check non-skel gens,
  // but currently we do not have a generic way to filter these out.
  const auto& y_cs = y.cons();
  const auto& x_gs = x.gens();
  for (const auto& c : y_cs) {
    if (not c.is_equality())
      continue;
    Index_Set c_nz;
    if (sp::is_sparse(c, c_nz)) {
      for (const auto& g : x_gs) {
        if (sp::violated_by(c_nz, c, g))
          return true;
      }
    } else {
      for (const auto& g : x_gs) {
        if (sp::violated_by(c, g))
          return true;
      }
    }
  }
  return false;
}

template <typename PH>
inline bool
widening_preamble(PH& x, const PH& y, Widen_Spec w_spec) {
  x.minimize();
  y.minimize();
  const bool safe_widen = (w_spec == Widen_Spec::SAFE);
  // Filter out empty and zero dim special cases.
  if (safe_widen) {
    // Safe widening spec: no containment precondition.
    if (y.is_empty())
      return true;
    if (x.is_empty()) {
      x = y;
      return true;
    }
    if (x.space_dim() == 0) {
      assert(x.is_universe() && y.is_universe());
      return true;
    }
  } else {
    // Risky widening spec: x contains y.
    assert(x.contains(y));
    if (x.space_dim() == 0 || x.is_empty() || y.is_empty())
      return true;
  }

  // Check if joining y with x would change (i.e., increase)
  // the affine dimension of x; this is not formally part of
  // a generic widening specification, but it is more precise
  // (and maybe also efficient).
  if (safe_widen) {
    if (increases_affine_dim(x, y)) {
      x.join_assign(y);
      return true;
    }
  } else {
    // Risky widen: can directly check affine dimension.
    if (x.affine_dim() > y.affine_dim())
      return true;
  }

  // Not a special case.
  return false;
}

inline Index_Set
valid_upto_cons(const Poly_Impl& ph, const Cons* upto_ptr) {
  assert(ph.is_minimized());
  Index_Set valid;
  if (not upto_ptr)
    return valid;
  const auto& upto_cs = *upto_ptr;
  for (auto i : bwd_index_range(upto_cs)) {
    if (sp::satisfied_by_all(ph.topol, upto_cs[i], ph.gs))
      valid.set(i);
  }
  return valid;
}

template <typename PH>
Index_Set
valid_upto_cons(const PH& ph, const Cons* upto_ptr) {
  Index_Set valid;
  if (not upto_ptr)
    return valid;
  const auto& upto_cs = *upto_ptr;
  for (auto i : bwd_index_range(upto_cs)) {
    auto rel = ph.relation_with(upto_cs[i]);
    if (rel.implies(Poly_Con_Rel::is_included()))
      valid.set(i);
  }
  return valid;
}

template <typename PH>
void
add_valid_upto_cons(PH& ph, const Index_Set& valid, const Cons* upto_ptr) {
  if (not upto_ptr)
    return;
  const auto& upto_cs = *upto_ptr;
  for (auto i : valid)
    ph.add_con(upto_cs[i]);
}

} // namespace detail

namespace h79_widen {

void h79_widen(Poly_Impl& x, const Poly_Impl& y, Widen_Spec w_spec);
void boxed_h79_widen(Poly_Impl& x, const Poly_Impl& y, Widen_Spec w_spec);

} // namespace h79_widen

namespace bhrz03_widen {

// The convergence certificate.
class Cert {
public:
  Cert() = default;
  Cert(const Cert&) = default;
  Cert& operator=(const Cert&) = default;
  Cert(Cert&&) noexcept = default;
  Cert& operator=(Cert&&) noexcept = default;
  ~Cert() = default;

  explicit Cert(const Poly& ph) : Cert(ph.impl()) {}
  explicit Cert(const Poly_Impl& ph);
  explicit Cert(const F_Poly& ph);
  template <typename PH>
  explicit Cert(const B_Wrap<PH>& ph) : Cert(ph.impl_poly()) {}

  int compare(const Cert& y) const;
  int compare(const Poly_Impl& ph) const;

  template <typename PH>
  int compare(const PH& ph) const { return compare(Cert(ph)); }

  // Returns true if ph is stabilizing wrt this certificate,
  // i.e., this is strictly bigger than ph in the ordering.
  template <typename PH>
  bool is_stabilizing(const PH& ph) const { return compare(ph) == 1; }

  bool check_inv() const;

private:
  // Note: default values are for the universe 0-dim polyhedron.
  dim_type affine_dim = 0;
  dim_type lin_space_dim = 0;
  // Non-redundant skeleton constraints.
  dim_type num_sk_cons = 0;
#ifndef NDEBUG
  // Just for debugging purposes.
  dim_type equalities = 0;
#endif
  /*
    A map providing a partial description of the supports of non-redundant
    strict inequalities; the meaning of a pair <j, k> is that there are
    k non-redundant strict inequalities whose support has cardinality j.
  */
  std::map<dim_type, dim_type> support_cards;
  // Non-redundant skeleton (closure) points.
  dim_type num_sk_points = 1;
  /*
    A vector containing, for each index `0 <= i < space_dim',
    the number of non-redundant rays in a generator system of the
    polyhedron having exactly `i' null coordinates.
  */
  std::vector<dim_type> num_rays_null_coord;
};

NOTHROW_MOVES(Cert);

void bhrz03_widen(Poly_Impl& x, const Poly_Impl& y, Widen_Spec w_spec);

} // namespace bhrz03_widen

} // namespace pplite

#endif // !defined(pplite_Poly_widen_hh)

