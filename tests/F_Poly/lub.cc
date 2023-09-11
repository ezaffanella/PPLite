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

#include "pplite_test.hh"

namespace pplite {
namespace detail {

using Blocks = F_Poly::Blocks;

Blocks least_upper_bound(const Blocks&, const Blocks&);

} // namespace detail
} // namespace pplite

bool
test01() {
  F_Poly::Blocks bs1 = {
    { 0 },
    { 3 },
    { 5 },
    { 10 },
    { 12 },
    { 6 },
    { 4 },
    { 1 },
    { 8, 9, 7 },
    { 11 },
    { 2 }
  };

  F_Poly::Blocks bs2 = {
    { 0, 8 },
    { 3 },
    { 5 },
    { 10 },
    { 11 },
    { 7, 6 },
    { 2 },
    { 1 },
    { 4 },
    { 12 },
    { 9 }
  };

  F_Poly::Blocks kr = {
    { 0, 8, 9, 7, 6 },
    { 3 },
    { 5 },
    { 10 },
    { 12 },
    { 4 },
    { 1 },
    { 11 },
    { 2 }
  };

  F_Poly::Blocks lub = detail::least_upper_bound(bs1, bs2);

  return lub == kr;
}

BEGIN_MAIN
  DO_TEST(test01);
END_MAIN
