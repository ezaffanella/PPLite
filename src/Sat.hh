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

#ifndef pplite_Sat_hh
#define pplite_Sat_hh 1

#include "globals.hh"
#include "Bits.hh"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <utility>
#include <vector>

namespace pplite {

class Sat {
public:
  using Row = Bits;
  using Rows = std::vector<Row>;

  struct Impl {
    Rows rows;
    dim_type row_size;
    Impl() noexcept : rows(), row_size(0) {}
    Impl(dim_type n_rows, dim_type n_cols) : rows(n_rows), row_size(n_cols) {}
  };

  Impl& impl() { return impl_; }
  const Impl& impl() const { return impl_; }

  Sat() = default;
  Sat(const Sat& y) = default;
  Sat(Sat&& y) = default;
  Sat& operator=(const Sat& y) = default;
  Sat& operator=(Sat&& y) = default;
  ~Sat() = default;

  Sat(dim_type n_rows, dim_type n_cols) : impl_(n_rows, n_cols) {}

  void swap(Sat& y) noexcept {
    using std::swap;
    swap(impl().rows, y.impl().rows);
    swap(impl().row_size, y.impl().row_size);
  }

  const Rows& rows() const { return impl().rows; }
  dim_type num_rows() const { return impl().rows.size(); }
  dim_type num_cols() const { return impl().row_size; }
  bool empty() const { return num_rows() == 0; }
  bool check_inv() const;

  Rows::const_iterator cbegin() const { return impl().rows.cbegin(); }
  Rows::const_iterator cend() const { return impl().rows.cend(); }

  Bits& operator[](dim_type k) {
    assert(k < num_rows());
    return impl().rows[k];
  }
  const Bits& operator[](dim_type k) const {
    assert(k < num_rows());
    return rows()[k];
  }

  void clear() { impl() = Impl(); }

  Sat transpose() const;

  void add_row(Row&& row) { impl().rows.push_back(std::move(row)); }

  void add_cols(dim_type n) { impl().row_size += n; }

  void remove_trailing_rows(dim_type n) {
    assert(n <= num_rows());
    impl().rows.resize(num_rows() - n);
  }
  void remove_trailing_columns(dim_type n) {
    assert(n <= num_cols());
    impl().row_size -= n;
    assert(check_inv());
  }
  void resize(dim_type new_nrows, dim_type new_ncols);

  void ascii_dump(std::ostream& s) const;
  bool ascii_load(std::istream& s);

  void print(std::ostream& os) const { ascii_dump(os); }
  void print() const { print(std::cout); }

private:
  Impl impl_;

}; // class Sat

NOTHROW_DEFAULT_AND_MOVES(Sat);

inline void swap(Sat& x, Sat& y) noexcept { x.swap(y); }

inline bool
operator==(const Sat& x, const Sat& y) {
  return x.num_rows() == y.num_rows()
    && x.num_cols() == y.num_cols()
    && x.rows() == y.rows();
}

inline bool
operator!=(const Sat& x, const Sat& y) { return !(x == y); }

} // namespace pplite

#endif // !defined(pplite_Sat_hh)
