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
#include "Con.hh"
#include "ascii_dump_load.hh"
#include <string>

namespace pplite {

bool
Con::is_tautological() const {
  switch (type()) {
  case EQUALITY:
    if (inhomo_term() != 0)
      return false;
    break;
  case NONSTRICT_INEQUALITY:
    if (inhomo_term() < 0)
      return false;
    break;
  case STRICT_INEQUALITY:
    if (inhomo_term() <= 0)
      return false;
    break;
  }
  return linear_expr().is_zero();
}

bool
Con::is_inconsistent() const {
  switch (type()) {
  case EQUALITY:
    if (inhomo_term() == 0)
      return false;
    break;
  case NONSTRICT_INEQUALITY:
    if (inhomo_term() >= 0)
      return false;
    break;
  case STRICT_INEQUALITY:
    if (inhomo_term() > 0)
      return false;
    break;
  }
  return linear_expr().is_zero();
}

bool
Con::check_inv() const {
  if (!check_strong_normalized()) {
#ifndef NDEBUG
    std::cerr << "Constraints should be strongly normalized!"
              << std::endl;
#endif
    return false;
  }
  return true;
}

void
Con::print(std::ostream& s) const {
  linear_expr().print(s);
  s << " ";
  dump_type(s);
  s << " ";
  switch (sgn(inhomo_term())) {
  case -1:
    (-inhomo_term()).print(s);
    break;
  case 0:
    s << "0";
    break;
  case 1:
    s << "-";
    inhomo_term().print(s);
    break;
  }
}

void
Con::dump_type(std::ostream& s) const {
  switch (type()) {
  case EQUALITY:
    s << "=";
    break;
  case NONSTRICT_INEQUALITY:
    s << ">=";
    break;
  case STRICT_INEQUALITY:
    s << ">";
    break;
  }
}

bool
Con::load_type(std::istream& s) {
  std::string str;
  if (!(s >> str))
    return false;
  if (str == "=")
    set_type(EQUALITY);
  else if (str == ">=")
    set_type(NONSTRICT_INEQUALITY);
  else if (str == ">")
    set_type(STRICT_INEQUALITY);
  else
    return false;
  return true;
}

void
Con::ascii_dump(std::ostream& s) const {
  dump_type(s);
  s << " : ";
  linear_expr().ascii_dump(s);
  s << " : ";
  inhomo_term().print(s);
  s << "\n";
}

bool
Con::ascii_load(std::istream& s) {
  return load_type(s)
    && ascii_load_string(s, ":")
    && impl().expr.ascii_load(s)
    && ascii_load_string(s, ":")
    && impl().inhomo.ascii_load(s);
}

namespace IO_Operators {

std::ostream&
operator<<(std::ostream& s, const Con& c) {
  bool first = true;
  for (auto i : dim_range(c)) {
    if (c.linear_expr().get(i) == 0)
      continue;
    Integer ci = c.linear_expr().get(i);
    if (first)
      first = false;
    else {
      if (ci > 0)
        s << " + ";
      else {
        s << " - ";
        neg_assign(ci);
      }
    }
    if (ci == -1)
      s << "-";
    else if (ci != 1)
      s << ci << "*";
    s << Var(i);
  }
  if (first)
    s << "0";
  const char* rel = nullptr;
  switch (c.type()) {
  case Con::EQUALITY:
    rel = " = ";
    break;
  case Con::NONSTRICT_INEQUALITY:
    rel = " >= ";
    break;
  case Con::STRICT_INEQUALITY:
    rel = " > ";
    break;
  }
  s << rel << -c.inhomo_term();
  return s;
}

std::ostream&
operator<<(std::ostream& s, const Con::Type& t) {
  const char* n = nullptr;
  switch (t) {
  case Con::EQUALITY:
    n = "EQUALITY";
    break;
  case Con::NONSTRICT_INEQUALITY:
    n = "NONSTRICT_INEQUALITY";
    break;
  case Con::STRICT_INEQUALITY:
    n = "STRICT_INEQUALITY";
    break;
  }
  s << n;
  return s;
}

std::ostream&
operator<<(std::ostream& s, const Cons& cs) {
  auto i = cs.begin(), i_end = cs.end();
  if (i == i_end) {
    s << "true";
    return s;
  }
  s << *i;
  for (++i; i != i_end; ++i)
    s << ", " << *i;
  return s;
}

} // namespace IO_Operators

} // namespace pplite

