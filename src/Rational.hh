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

#ifndef pplite_Rational_hh
#define pplite_Rational_hh 1

#include "globals.hh"
#include "Integer.hh"
#include "Rational_fwd.hh"

#if PPLITE_USE_FLINT_INTEGERS

#include "FLINT_Rational.hh"

#else // !PPLITE_USE_FLINT_INTEGERS

#include "GMP_Rational.hh"

#endif // !PPLITE_USE_FLINT_INTEGERS

#include <vector>

namespace pplite {

inline std::ostream&
operator<<(std::ostream& s, const Rational& r) {
  r.print(s);
  return s;
}

inline std::istream&
operator>>(std::istream& s, Rational& r) {
  if (!r.ascii_load(s))
    s.setstate(std::ios::failbit);
  return s;
}

using Rationals = std::vector<Rational>;

inline void
rationals_to_integers(const Rationals& src,
                      Integers& dst, Integer& denom) {
  assert(src.size() == dst.size());
  denom = lcm_dens(src);
  for (auto i : index_range(src))
    exact_div_assign(dst[i], denom, src[i]);
}

} // namespace pplite

#endif // !defined(pplite_Rational_hh)
