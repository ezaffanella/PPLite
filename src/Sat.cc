/* PPLite: a lightweight library for convex polyhedra derived from PPL.
   Copyright (C) 2001-2010 Roberto Bagnara <bagnara@cs.unipr.it>
   Copyright (C) 2010-2018 BUGSENG srl (http://bugseng.com)
   Copyright (C) 2018-2024 Enea Zaffanella <enea.zaffanella@unipr.it>

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
#include "Sat.hh"

#include <string>

namespace pplite {

void
Sat::resize(dim_type new_nrows, dim_type new_ncols) {
  Impl& x = impl();
  x.rows.resize(new_nrows);
  if (new_ncols < x.row_size) {
    for (auto& r : x.rows)
      r.reset_from(new_ncols);
  }
  x.row_size = new_ncols;
  assert(check_inv());
}

Sat
Sat::transpose() const {
  const auto& x = *this;
  Sat res(x.num_cols(), x.num_rows());
  for (auto i : bwd_range(x.num_rows()))
    // for (auto j = x[i].last(); j != Row::end_pos; j = x[i].prev(j))
    for (auto j = x[i].first(); j != Row::end_pos; j = x[i].next(j))
      res[j].set(i);
  assert(res.check_inv());
  return res;
}

void
Sat::ascii_dump(std::ostream& s) const {
  const char sep = ' ';
  s << num_rows() << sep << 'x' << sep << num_cols() << "\n";
  for (const auto& row : rows()) {
    for (auto j : range(num_cols())) {
      s << row[j] << sep;
    }
    s << "\n";
  }
}

bool
Sat::ascii_load(std::istream& s) {
  auto& x = *this;
  dim_type nrows;
  dim_type ncols;
  std::string str;
  if (!(s >> nrows))
    return false;
  if (!(s >> str) || str != "x")
    return false;
  if (!(s >> ncols))
    return false;
  resize(nrows, ncols);

  for (auto i : range(nrows)) {
    for (auto j : range(ncols)) {
      int bit;
      if (!(s >> bit))
        return false;
      if (bit != 0)
        x[i].set(j);
      else
        x[i].reset(j);
    }
  }
  assert(check_inv());
  return true;
}

bool
Sat::check_inv() const {
#ifndef NDEBUG
  using std::endl;
  using std::cerr;
#endif

  const auto& x = *this;
  for (auto i : bwd_range(num_rows())) {
    const auto& row = x[i];
    if (!row.check_inv()) {
#ifndef NDEBUG
      cerr << "sat[" << i << "] is broken" << endl;
#endif
      return false;
    }
    if (row.last() != Row::end_pos && row.last() >= num_cols()) {
#ifndef NDEBUG
      cerr << "Sat[" << i << "] is a row with too many bits!"
           << endl
           << "(row_size == " << impl().row_size
           << ", row.last() == " << row.last() << ")"
           << endl;
#endif
      return false;
    }
  }
  return true;
}

} // namespace pplite

