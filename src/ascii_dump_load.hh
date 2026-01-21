/* PPLite: a lightweight library for convex polyhedra derived from PPL.
   Copyright (C) 2001-2010 Roberto Bagnara <bagnara@cs.unipr.it>
   Copyright (C) 2010-2018 BUGSENG srl (http://bugseng.com)
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

#ifndef pplite_ascii_dump_load_hh
#define pplite_ascii_dump_load_hh 1

#include <cctype>
#include <gmp.h>
#include <iterator>
#include <iostream>
#include <string>
#include <vector>

namespace pplite {

inline void
ascii_load_skip_spaces(std::istream& s) {
  char ch;
  while (s.get(ch)) {
    if (not std::isspace(ch)) {
      s.unget();
      return;
    }
  }
}

inline const char*
mpz_to_string(const mpz_t mp) {
  // Create buffer (done only once).
  static PPLITE_TLS std::vector<char> buffer;
  // Resize buffer as needed (+2 for sign and terminator)
  const int base10 = 10;
  const size_t size = 2 + mpz_sizeinbase(mp, base10);
  buffer.resize(size);
  // Fill buffer with char* representation
  mpz_get_str(buffer.data(), base10, mp);
  return buffer.data();
}

inline const char*
read_mpz_as_string(std::istream& is) {
  // Create buffer (done only once).
  static PPLITE_TLS std::vector<char> buffer;
  buffer.resize(0);
  char ch;
  if (not is.get(ch))
    return nullptr;
  // Check for (optional) minus sign
  if (ch == '-') {
    buffer.push_back(ch);
    if (not is.get(ch)) {
      is.unget();
      return nullptr;
    }
  }
  // Check for (required) first digit
  if (std::isdigit(ch))
    buffer.push_back(ch);
  else {
    // Not a digit: error
    is.unget();
    return nullptr;
  }
  // Handle any other digit
  while (is.get(ch)) {
    if (std::isdigit(ch))
      buffer.push_back(ch);
    else {
      // Not a digit: end of input
      is.unget();
      break;
    }
  }
  // Push back the terminator
  buffer.push_back('\0');
  return buffer.data();
}

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
