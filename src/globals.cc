/* PPLite: a lightweight library for convex polyhedra derived from PPL.
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
#include "globals.hh"

namespace pplite::detail {

// Initialize (thread-local) globals.
PPLITE_TLS Topol default_topol = Topol::CLOSED;
PPLITE_TLS Widen_Spec widen_spec = Widen_Spec::SAFE;
PPLITE_TLS Widen_Impl widen_impl = Widen_Impl::H79;

} // namespace pplite::detail
