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
Con::print(std::ostream& os) const {
  linear_expr().print(os);
  os << " ";
  dump_type(os);
  os << " ";
  switch (sgn(inhomo_term())) {
  case -1:
    (-inhomo_term()).print(os);
    break;
  case 0:
    os << "0";
    break;
  case 1:
    os << "-";
    inhomo_term().print(os);
    break;
  }
}

void
Con::dump_type(std::ostream& os) const {
  switch (type()) {
  case EQUALITY:
    os << "=";
    break;
  case NONSTRICT_INEQUALITY:
    os << ">=";
    break;
  case STRICT_INEQUALITY:
    os << ">";
    break;
  }
}

bool
Con::load_type(std::istream& is) {
  std::string str;
  if (!(is >> str))
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
Con::ascii_dump(std::ostream& os) const {
  dump_type(os);
  os << " : ";
  linear_expr().ascii_dump(os);
  os << " : ";
  inhomo_term().print(os);
  os << "\n";
}

bool
Con::ascii_load(std::istream& is) {
  return load_type(is)
    && ascii_load_string(is, ":")
    && impl().expr.ascii_load(is)
    && ascii_load_string(is, ":")
    && impl().inhomo.ascii_load(is);
}

namespace IO_Operators {

std::ostream&
operator<<(std::ostream& os, const Con& c) {
  bool first = true;
  for (auto i : dim_range(c)) {
    if (c.linear_expr().get(i) == 0)
      continue;
    Integer ci = c.linear_expr().get(i);
    if (first)
      first = false;
    else {
      if (ci > 0)
        os << " + ";
      else {
        os << " - ";
        neg_assign(ci);
      }
    }
    if (ci == -1)
      os << "-";
    else if (ci != 1)
      os << ci << "*";
    os << Var(i);
  }
  if (first)
    os << "0";
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
  os << rel << -c.inhomo_term();
  return os;
}

std::ostream&
operator<<(std::ostream& os, const Con::Type& t) {
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
  os << n;
  return os;
}

std::ostream&
operator<<(std::ostream& os, const Cons& cs) {
  auto i = cs.begin(), i_end = cs.end();
  if (i == i_end) {
    os << "true";
    return os;
  }
  os << *i;
  for (++i; i != i_end; ++i)
    os << ", " << *i;
  return os;
}

} // namespace IO_Operators

} // namespace pplite

