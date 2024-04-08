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

#ifndef pplite_PolySet_templ_hh
#define pplite_PolySet_templ_hh

#include "PolySet.hh"
#include <algorithm>
#include <iterator>
#include <type_traits>

namespace pplite {

// Checking if a type is instantiated from PolySet.
template <typename>
struct is_PolySet : public std::false_type { };
template <typename PH>
struct is_PolySet<PolySet<PH>> : public std::true_type { };

template <typename PH>
template <typename Iter>
void
PolySet<PH>::add_cons(Iter first, Iter last) {
  for (auto& d : seq())
    d.add_cons(first, last);
  clear_reduced();
  assert(check_inv());
}

template <typename PH>
template <typename Iter>
void
PolySet<PH>::add_gens(Iter first, Iter last) {
  if (seq().empty())
    seq().emplace_back(space_dim(), topology(), Spec_Elem::EMPTY);
  for (auto& d : seq())
    d.add_gens(first, last);
  clear_reduced();
  assert(check_inv());
}

template <typename PH>
template <typename Iter>
void
PolySet<PH>::unconstrain(Iter first, Iter last) {
  for (auto& d : seq())
    d.unconstrain(first, last);
  clear_reduced();
  assert(check_inv());
}

template <typename PH>
template <typename Iter>
void
PolySet<PH>::remove_space_dims(Iter first, Iter last) {
  auto num_removed = std::distance(first, last);
  if (num_removed == 0)
    return;
  for (auto& d : seq()) {
    d.remove_space_dims(first, last);
    clear_reduced();
  }
  impl().dim -= num_removed;
  assert(check_inv());
}

} // namespace pplite

#endif //!defined(pplite_PolySet_templ_hh)
