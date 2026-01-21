/* PPLite: a lightweight library for convex polyhedra derived from PPL.
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
#include "Abs_Poly_Adapter.hh"

namespace pplite {
namespace dynamic {

// Definitions of explicit instantiations of the available concrete classes.
template class Abs_Poly_Adapter<Poly>;
template class Abs_Poly_Adapter<Stats<Poly>>;
template class Abs_Poly_Adapter<B_Poly>;
template class Abs_Poly_Adapter<Stats<B_Poly>>;
template class Abs_Poly_Adapter<F_Poly>;
template class Abs_Poly_Adapter<Stats<F_Poly>>;
template class Abs_Poly_Adapter<U_Poly>;
template class Abs_Poly_Adapter<Stats<U_Poly>>;
template class Abs_Poly_Adapter<UF_Poly>;
template class Abs_Poly_Adapter<Stats<UF_Poly>>;

template class Abs_Poly_Adapter<P_Set>;
template class Abs_Poly_Adapter<Stats<P_Set>>;
template class Abs_Poly_Adapter<FP_Set>;
template class Abs_Poly_Adapter<Stats<FP_Set>>;

PPLITE_TLS Abs_Poly::Kind default_poly_kind = Abs_Poly::Kind::POLY;

bool
abs_poly_name_to_kind(const char* cname,
                      Abs_Poly::Kind& kind, bool& requires_stats) {
  std::string sname = cname;
  for (auto i : range(abs_poly_kind_size)) {
    if (sname == abs_poly_kind_names[i]) {
      kind = static_cast<Abs_Poly::Kind>(i);
      requires_stats = (sname.find("Stats") != std::string::npos);
      return true;
    }
  }
  return false;
}

bool
set_default_poly_kind(const char* cname, bool noisy) {
  auto kind = Abs_Poly::Kind::POLY;
  bool requires_stats = false;
  if (abs_poly_name_to_kind(cname, kind, requires_stats)) {
    default_poly_kind = kind;
    if (requires_stats)
      set_noisy_stats(true);
    if (noisy) {
      std::string msg = "Abs_Poly kind set to ";
      msg += cname;
      std::string marker = std::string(msg.size(), '=');
      std::cout << marker
                << "\n" << msg << "\n"
                << marker
                << std::endl;
    }
    return true;
  }
  if (noisy) {
    std::string msg = "Unknown Abs_Poly kind name '";
    msg += cname;
    msg += "'.\nValid values:";
    bool comma = false;
    for (const auto& name : abs_poly_kind_names) {
      msg += (comma ? ", " : " ") + std::string(name);
      comma = true;
    }
    std::cout << msg << std::endl;
  }
  return false;
}

} // namespace dynamic
} // namespace pplite
