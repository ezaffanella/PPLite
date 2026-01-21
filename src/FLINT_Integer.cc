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

#include "FLINT_Integer.hh"
#include "ascii_dump_load.hh"

#include <cctype>
#include <string>

namespace pplite {

void
FLINT_Integer::print(std::ostream& os, const fmpz_t mp) {
  const auto& d = *mp;
  if (COEFF_IS_MPZ(d)) {
    // Get the mpz_t.
    auto val = COEFF_TO_PTR(d);
    os << mpz_to_string(val);
  } else {
    // Here the value fits a signed integer
    assert(fmpz_fits_si(mp));
    os << fmpz_get_si(mp);
  }
}

bool
FLINT_Integer::read(std::istream& is, fmpz_t mp) {
  const char* buffer = read_mpz_as_string(is);
  if (buffer == nullptr)
    return false;
  int res = fmpz_set_str(mp, buffer, 10);
  return (res == 0);
}

void
FLINT_Integer::ascii_dump(std::ostream& s) const {
  print(s, mp);
}

bool
FLINT_Integer::ascii_load(std::istream& is) {
  ascii_load_skip_spaces(is);
  return read(is, mp);
}

} // namespace pplite

