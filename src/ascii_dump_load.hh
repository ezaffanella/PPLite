/* PPLite: a lightweight library for convex polyhedra derived from PPL.
   Copyright (C) 2001-2010 Roberto Bagnara <bagnara@cs.unipr.it>
   Copyright (C) 2010-2018 BUGSENG srl (http://bugseng.com)
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

#ifndef pplite_ascii_dump_load_hh
#define pplite_ascii_dump_load_hh 1

#include <string>
#include <iterator>
#include <iostream>

namespace pplite {

inline bool
ascii_load_string(std::istream& s, const std::string& str) {
  std::string loaded;
  return (s >> loaded) && (loaded == str);
}

template <typename Iter>
inline void
ascii_dump(std::ostream& s, Iter first, Iter last) {
  for ( ; first != last; ++first)
    first->ascii_dump(s);
}

template <typename Rows>
inline void
ascii_dump_all(std::ostream& s, const Rows& rows) {
  ascii_dump(s, std::begin(rows), std::end(rows));
}

} // namespace pplite

#endif // !defined(pplite_ascii_dump_load_hh)
