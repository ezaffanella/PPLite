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
Gen::dump_type(std::ostream& s) const {
  switch (type()) {
  case LINE:
    s << "L";
    break;
  case RAY:
    s << "R";
    break;
  case POINT:
    s << "P";
    break;
  case CLOSURE_POINT:
    s << "C";
    break;
  }
}

bool
Gen::load_type(std::istream& s) {
  std::string str;
  if (!(s >> str))
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
Gen::ascii_dump(std::ostream& s) const {
  dump_type(s);
  s << " : ";
  linear_expr().ascii_dump(s);
  s << " : ";
  impl().inhomo.print(s);
  s << "\n";
}

bool
Gen::ascii_load(std::istream& s) {
  return load_type(s)
    && ascii_load_string(s, ":")
    && impl().expr.ascii_load(s)
    && ascii_load_string(s, ":")
    && impl().inhomo.ascii_load(s);
}

void
Gen::print(std::ostream& s) const {
  bool need_div = false;
  switch (type()) {
  case LINE:
    s << "l(";
    break;
  case RAY:
    s << "r(";
    break;
  case POINT:
    s << "p(";
    need_div = (divisor() != 1);
    break;
  case CLOSURE_POINT:
    s << "c(";
    need_div = (divisor() != 1);
    break;
  }
  if (need_div)
    s << "(";
  linear_expr().print(s);
  if (need_div)
    s << ")/" << divisor();
  s << ")";
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
operator<<(std::ostream& s, const Gen::Type& t) {
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
  s << n;
  return s;
}

std::ostream&
operator<<(std::ostream& s, const Gens& gs) {
  auto i = gs.begin(), i_end = gs.end();
  if (i == i_end) {
    s << "false";
    return s;
  }
  s << *i;
  for (++i; i != i_end; ++i)
    s << ", " << *i;
  return s;
}

} // namespace IO_Operators

} // namespace pplite

