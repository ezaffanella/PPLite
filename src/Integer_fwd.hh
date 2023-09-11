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

#ifndef pplite_Integer_fwd_hh
#define pplite_Integer_fwd_hh 1

namespace pplite {

#if PPLITE_USE_FLINT_INTEGERS

class FLINT_Integer;
using Integer = FLINT_Integer;

#else // !PPLITE_USE_FLINT_INTEGER

class GMP_Integer;
using Integer = GMP_Integer;

#endif // !PPLITE_USE_FLINT_INTEGERS

} // namespace pplite

#endif // !defined(pplite_Integer_fwd_hh)
