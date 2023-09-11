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

#include "FLINT_Integer.hh"
#include <gmpxx.h>

namespace pplite {

void
FLINT_Integer::ascii_dump(std::ostream& s) const {
  print(s);
}

bool
FLINT_Integer::ascii_load(std::istream& is) {
  // flint has no operator>> : use the one in gmpxx.
  mpz_class z;
  if (!(is >> z))
    return false;
  fmpz_set_mpz(impl(), z.get_mpz_t());
  return true;
}

} // namespace pplite

