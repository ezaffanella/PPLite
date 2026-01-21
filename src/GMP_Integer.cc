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

#include "GMP_Integer.hh"
#include "ascii_dump_load.hh"

namespace pplite {

void
GMP_Integer::print(std::ostream& os, const mpz_t mp) {
  os << mpz_to_string(mp);
}

bool
GMP_Integer::read(std::istream& is, mpz_t mp) {
  const char* buffer = read_mpz_as_string(is);
  if (buffer == nullptr)
    return false;
  int res = mpz_set_str(mp, buffer, 10);
  return (res == 0);
}

void
GMP_Integer::ascii_dump(std::ostream& s) const {
  print(s, mp);
}

bool
GMP_Integer::ascii_load(std::istream& is) {
  ascii_load_skip_spaces(is);
  return read(is, mp);
}

} // namespace pplite

