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

#ifndef pplite_Integer_hh
#define pplite_Integer_hh 1

#include "Integer_fwd.hh"

#include "globals.hh"
#include "ascii_dump_load.hh"
#include <vector>

#if PPLITE_USE_FLINT_INTEGERS

#include "FLINT_Integer.hh"

#else // !PPLITE_USE_FLINT_INTEGERS

#include "GMP_Integer.hh"

#endif // !PPLITE_USE_FLINT_INTEGERS

namespace pplite {

inline bool using_FLINT() {
  return PPLITE_USE_FLINT_INTEGERS;
}

inline void
get_coprimes(const Integer& x, const Integer& y,
             Integer& co_x, Integer& co_y) {
  // Note: use co_y to store gcd to avoid temporary;
  // exact_div_assign allows for argument aliasing.
  Integer& gcd = co_y;
  gcd_assign(gcd, x, y);
  exact_div_assign(co_x, x, gcd);
  exact_div_assign(co_y, y, gcd);
}

inline std::ostream&
operator<<(std::ostream& s, const Integer& i) {
  i.print(s);
  return s;
}

inline std::istream&
operator>>(std::istream& s, Integer& i) {
  if (s.flags() & std::ios::skipws)
    ascii_load_skip_spaces(s);
  if (not i.read(s))
    s.setstate(std::ios::failbit);
  return s;
}

using Integers = std::vector<Integer>;

} // namespace pplite

#endif // !defined(pplite_Integer_hh)
