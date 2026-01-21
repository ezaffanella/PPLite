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

#ifndef pplite_Poly_Rel_hh
#define pplite_Poly_Rel_hh 1

#include "globals.hh"
#include <iostream>

namespace pplite {

class Poly_Con_Rel {
private:
  using Impl = unsigned int;

  static constexpr Impl NOTHING             = 0U;
  static constexpr Impl IS_DISJOINT         = 1U << 0;
  static constexpr Impl STRICTLY_INTERSECTS = 1U << 1;
  static constexpr Impl IS_INCLUDED         = 1U << 2;
  static constexpr Impl SATURATES           = 1U << 3;

  static constexpr Impl EVERYTHING
  = IS_DISJOINT
  | STRICTLY_INTERSECTS
  | IS_INCLUDED
  | SATURATES;

  Impl flags;

  // Not explicit: behaves as implicit conversion.
  constexpr Poly_Con_Rel(Impl mask) : flags(mask) {}
  static bool implies(Impl x, Impl y) { return (x & y) == y; }

public:
  constexpr Poly_Con_Rel() noexcept : flags(NOTHING) {}

  Impl& impl() { return flags; }
  Impl impl() const { return flags; }

  static constexpr Poly_Con_Rel nothing() { return NOTHING; }
  static constexpr Poly_Con_Rel is_disjoint() { return IS_DISJOINT; }
  static constexpr Poly_Con_Rel
  strictly_intersects() { return STRICTLY_INTERSECTS; }
  static constexpr Poly_Con_Rel is_included() { return IS_INCLUDED; }
  static constexpr Poly_Con_Rel saturates() { return SATURATES; }

  bool implies(const Poly_Con_Rel& y) const {
    return implies(flags, y.flags);
  }

  void ascii_dump(std::ostream& s) const;
  void print(std::ostream& os) const { ascii_dump(os); }
  void print() const { print(std::cout); }

}; // class Poly_Con_Rel

class Poly_Gen_Rel {
private:
  using Impl = unsigned int;

  static constexpr Impl NOTHING  = 0U;
  static constexpr Impl SUBSUMES = 1U << 0;

  static constexpr Impl EVERYTHING = SUBSUMES;

  Impl flags;

  // Not explicit: behaves as implicit conversion.
  constexpr Poly_Gen_Rel(Impl mask) : flags(mask) {}
  static bool implies(Impl x, Impl y) { return (x & y) == y; }

public:
  constexpr Poly_Gen_Rel() noexcept : flags(NOTHING) {}

  Impl& impl() { return flags; }
  Impl impl() const { return flags; }

  static constexpr Poly_Gen_Rel nothing() { return NOTHING; }
  static constexpr Poly_Gen_Rel subsumes() { return SUBSUMES; }

  bool implies(const Poly_Gen_Rel& y) const {
    return implies(impl(), y.impl());
  }

  void ascii_dump(std::ostream& s) const;
  void print(std::ostream& os) const { ascii_dump(os); }
  void print() const { print(std::cout); }

}; // class Poly_Gen_Rel

NOTHROW_DEFAULT_AND_MOVES(Poly_Con_Rel);
NOTHROW_DEFAULT_AND_MOVES(Poly_Gen_Rel);

inline bool
operator==(const Poly_Con_Rel& x, const Poly_Con_Rel& y) {
  return x.impl() == y.impl();
}
inline bool
operator==(const Poly_Gen_Rel& x, const Poly_Gen_Rel& y) {
  return x.impl() == y.impl();
}

inline bool
operator!=(const Poly_Con_Rel& x, const Poly_Con_Rel& y) {
  return !(x == y);
}
inline bool
operator!=(const Poly_Gen_Rel& x, const Poly_Gen_Rel& y) {
  return !(x == y);
}

inline Poly_Con_Rel
operator&&(Poly_Con_Rel x, Poly_Con_Rel y) {
  x.impl() |= y.impl();
  return x;
}
inline Poly_Gen_Rel
operator&&(Poly_Gen_Rel x, Poly_Gen_Rel y) {
  x.impl() |= y.impl();
  return x;
}

inline Poly_Con_Rel
operator-(Poly_Con_Rel x, Poly_Con_Rel y) {
  x.impl() &= ~y.impl();
  return x;
}
inline Poly_Gen_Rel
operator-(Poly_Gen_Rel x, Poly_Gen_Rel y) {
  x.impl() &= ~y.impl();
  return x;
}

namespace IO_Operators {

inline std::ostream&
operator<<(std::ostream& s, const Poly_Con_Rel& r) {
  r.ascii_dump(s);
  return s;
}
inline std::ostream&
operator<<(std::ostream& s, const Poly_Gen_Rel& r) {
  r.ascii_dump(s);
  return s;
}

} // namespace IO_Operators

} // namespace pplite

#endif // !defined(pplite_Poly_Rel_hh)
