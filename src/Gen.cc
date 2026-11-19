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

#include "pplite-config.h"
#include "Gen.hh"
#include "ascii_dump_load.hh"
#include <string>

namespace pplite {

bool
Gen::check_inv() const {
  if (!check_strong_normalized()) {
#ifndef NDEBUG
    std::cerr << "Generators should be strongly normalized!"
              << std::endl;
#endif
    return false;
  }

  if (is_line_or_ray()) {
    if (impl().inhomo != 0) {
#ifndef NDEBUG
      std::cerr << "Lines and rays must have a zero divisor!"
                << std::endl;
#endif
      return false;
    }
    if (linear_expr().is_zero()) {
#ifndef NDEBUG
      std::cerr << "The origin of the vector space cannot be a line or a ray!"
                << std::endl;
#endif
      return false;
    }
  }
  else {
    if (divisor() <= 0) {
#ifndef NDEBUG
      std::cerr << "Points and closure points must have a positive divisor!"
                << std::endl;
#endif
      return false;
    }
  }

  // All tests passed.
  return true;
}

void
Gen::dump_type(std::ostream& os) const {
  switch (type()) {
  case LINE:
    os << "L";
    break;
  case RAY:
    os << "R";
    break;
  case POINT:
    os << "P";
    break;
  case CLOSURE_POINT:
    os << "C";
    break;
  }
}

bool
Gen::load_type(std::istream& is) {
  std::string str;
  if (!(is >> str))
    return false;
  if (str == "P")
    set_type(POINT);
  else if (str == "R")
    set_type(RAY);
  else if (str == "L")
    set_type(LINE);
  else if (str == "C")
    set_type(CLOSURE_POINT);
  else
    return false;
  return true;
}

void
Gen::ascii_dump(std::ostream& os) const {
  dump_type(os);
  os << " : ";
  linear_expr().ascii_dump(os);
  os << " : ";
  impl().inhomo.print(os);
  os << "\n";
}

bool
Gen::ascii_load(std::istream& is) {
  return load_type(is)
    && ascii_load_string(is, ":")
    && impl().expr.ascii_load(is)
    && ascii_load_string(is, ":")
    && impl().inhomo.ascii_load(is);
}

void
Gen::print(std::ostream& os) const {
  bool need_div = false;
  switch (type()) {
  case LINE:
    os << "l(";
    break;
  case RAY:
    os << "r(";
    break;
  case POINT:
    os << "p(";
    need_div = (divisor() != 1);
    break;
  case CLOSURE_POINT:
    os << "c(";
    need_div = (divisor() != 1);
    break;
  }
  if (need_div)
    os << "(";
  linear_expr().print(os);
  if (need_div)
    os << ")/" << divisor();
  os << ")";
}

namespace detail {

void
add_as_rays(Gens gens, Gens& rays) {
  for (auto& g : gens) {
    if (g.is_point() || g.is_closure_point()) {
      // Skip the origin (not a valid ray).
      if (g.linear_expr().is_zero())
        continue;
      // Make it a ray.
      g.impl().inhomo = 0;
      g.impl().expr.normalize(g.impl().inhomo);
      g.impl().type = Gen::RAY;
    }
    rays.push_back(std::move(g));
  }
}

} // namespace detail

namespace IO_Operators {

std::ostream&
operator<<(std::ostream& os, const Gen::Type& t) {
  const char* n = nullptr;
  switch (t) {
  case Gen::POINT:
    n = "POINT";
    break;
  case Gen::RAY:
    n = "RAY";
    break;
  case Gen::LINE:
    n = "LINE";
    break;
  case Gen::CLOSURE_POINT:
    n = "CLOSURE_POINT";
    break;
  }
  os << n;
  return os;
}

std::ostream&
operator<<(std::ostream& os, const Gens& gs) {
  auto i = gs.begin(), i_end = gs.end();
  if (i == i_end) {
    os << "false";
    return os;
  }
  os << *i;
  for (++i; i != i_end; ++i)
    os << ", " << *i;
  return os;
}

} // namespace IO_Operators

} // namespace pplite

