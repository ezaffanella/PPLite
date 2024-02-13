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

#include "FLINT_Rational.hh"
#include "ascii_dump_load.hh"

namespace pplite {

void
FLINT_Rational::print(std::ostream& os, const fmpq_t mp) {
  FLINT_Integer::print(os, fmpq_numref(mp));
  if (fmpz_is_one(fmpq_denref(mp)))
    return;
  os << "/";
  FLINT_Integer::print(os, fmpq_denref(mp));
}

bool
FLINT_Rational::read(std::istream& is, fmpq_t mp) {
  FLINT_Integer num;
  FLINT_Integer den = FLINT_Integer::one();
  if (not num.read(is))
    return false;
  // Note: the "/den" part is optional
  char ch;
  if (is.get(ch)) {
    if (ch == '/') {
      if (not den.read(is))
        return false;
    } else {
      is.unget();
    }
  }
  fmpq_set_fmpz_frac(mp, num.impl(), den.impl());
  return true;
}

void
FLINT_Rational::ascii_dump(std::ostream& s) const {
  print(s, mp);
}

bool
FLINT_Rational::ascii_load(std::istream& is) {
  ascii_load_skip_spaces(is);
  return read(is, mp);
}

} // namespace pplite

