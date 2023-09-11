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
#include "Var.hh"

#include <iostream>

namespace pplite {

void
Var_default_output_function(std::ostream& s, const Var v) {
  const dim_type varid = v.id();
  static constexpr char var_name_letters[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  constexpr dim_type num_letters = sizeof(var_name_letters) - 1;
  s << var_name_letters[varid % num_letters];
  if (const dim_type i = varid / num_letters) {
    s << i;
  }
}

PPLITE_TLS Output_Function<Var>
Var::output_function(Var_default_output_function);

} // namespace pplite
