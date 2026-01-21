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

#ifndef pplite_Output_Function_hh
#define pplite_Output_Function_hh 1

#include <functional>

template <typename T>
class Output_Function {
public:
  using ftype = void (std::ostream& s, const T& t);

  template <typename F>
  explicit Output_Function(F default_f) : current(default_f) {}

  template <typename F>
  void set(F f) { current = f; }

  std::function<ftype> get_current() const { return current; }

private:
  std::function<ftype> current;
}; // class Output_Function

#endif // !defined(pplite_Output_Function_hh)
