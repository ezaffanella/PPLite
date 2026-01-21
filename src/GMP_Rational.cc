/* PPLite: a lightweight library for convex polyhedra derived from PPL.
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

#include "Integer.hh"
#include "GMP_Rational.hh"
#include "ascii_dump_load.hh"

namespace pplite {

void
GMP_Rational::print(std::ostream& os, const mpq_t mp) {
  GMP_Integer::print(os, mpq_numref(mp));
  if (0 == mpz_cmp_si(mpq_denref(mp), 1))
    return;
  os << "/";
  GMP_Integer::print(os, mpq_denref(mp));
}

bool
GMP_Rational::read(std::istream& is, mpq_t mp) {
  GMP_Integer num;
  GMP_Integer den = GMP_Integer::one();
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
  mpq_set_num(mp, num.impl());
  mpq_set_den(mp, den.impl());
  mpq_canonicalize(mp);
  return true;
}

void
GMP_Rational::ascii_dump(std::ostream& s) const {
  print(s, mp);
}

bool
GMP_Rational::ascii_load(std::istream& is) {
  ascii_load_skip_spaces(is);
  return read(is, mp);
}

} // namespace pplite

