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

#ifndef pplite_Var_hh
#define pplite_Var_hh 1

#include "globals.hh"
#include "Bits.hh"
#include "Output_Function.hh"

#include <cassert>
#include <iosfwd>
#include <utility>
#include <vector>

namespace pplite {

class Var {
private:
  dim_type varid;

public:
  explicit Var(dim_type i) : varid(i) { assert(check_inv()); }

  Var(const Var&) = default;
  Var(Var&&) = default;
  Var& operator=(const Var&) = default;
  Var& operator=(Var&&) = default;
  ~Var() = default;

  void m_swap(Var& v) noexcept {
    using std::swap;
    swap(varid, v.varid);
  }

  dim_type id() const { return varid; }
  dim_type space_dim() const { return varid + 1; }
  bool check_inv() const { return varid >= 0; }

  static PPLITE_TLS Output_Function<Var> output_function;

  void print(std::ostream& os) const {
    output_function.get_current()(os, *this);
  }

}; // class Var

inline bool
less_than(Var v, Var w) { return v.id() < w.id(); }

inline void
swap(Var& v, Var& w) noexcept { v.m_swap(w); }

namespace IO_Operators {

inline std::ostream&
operator<<(std::ostream& os, const Var v) {
  v.print(os);
  return os;
}

} // namespace IO_Operators

NOTHROW_MOVES(Var);

using Vars = std::vector<Var>;

struct Vars_Set : public Bits {
  using Bits::Bits;
  void insert(Var var) { set(var.id()); }
};

} // namespace pplite

#endif // !defined(pplite_Var_hh)
