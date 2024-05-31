/* PPLite: a lightweight library for convex polyhedra derived from PPL.
   Copyright (C) 2018-2024 Enea Zaffanella <enea.zaffanella@unipr.it>

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

#ifndef pplite_globals_hh
#define pplite_globals_hh 1

#ifndef PPLITE_CONVERSION_TIME_STATS
// Produces number of calls and computation time for conversion:
// output occurs at end of program execution.
#define PPLITE_CONVERSION_TIME_STATS 0
#endif

#ifndef PPLITE_LOW_LEVEL_COUNTERS
// Enables collection of number of calls to low level operations
// (scalar products, linear combinations, saturation row inclusion tests,
// unions and population counts). See Low_Level_Stats.hh for hints on
// using LLOp_Stats::reset_stats() and LLOp_Stats::print_stats(std::ostream&).
#define PPLITE_LOW_LEVEL_COUNTERS 0
#endif

#ifndef PPLITE_NOISY_CONVERSION
// Produces noisy stats after each iteration of the conversion procedure.
#define PPLITE_NOISY_CONVERSION 0
#endif

#ifdef PPLITE_THREAD_SAFE
#define PPLITE_TLS thread_local
#else
#define PPLITE_TLS
#endif

#include <cstddef>
#include <limits>
#include <type_traits>
#include <vector>

#define PPLITE_UNIMPL                                             \
  do {                                                            \
    std::cerr << "\n\npplite: calling unimplemented function\n";  \
    std::cerr << "  " << __PRETTY_FUNCTION__ << "\n";             \
    assert(false);                                                \
    abort();                                                      \
  } while (false)

#define PPLITE_UNREACH                                            \
  do {                                                            \
    std::cerr << "\n\npplite: executing unreachable code\n";      \
    std::cerr << "  " << __PRETTY_FUNCTION__ << "\n";             \
    assert(false);                                                \
    abort();                                                      \
  } while (false)


#define NOTHROW_DEFAULT_CTOR(Type)                                  \
  static_assert(std::is_nothrow_default_constructible<Type>::value, \
                #Type " should be noexcept default constructible")

#define NOTHROW_MOVES(Type)                                      \
  static_assert(std::is_nothrow_move_constructible<Type>::value, \
                #Type " should be noexcept move constructible"); \
  static_assert(std::is_nothrow_move_assignable<Type>::value,    \
                #Type " should be noexcept move assignable")

#define NOTHROW_DEFAULT_AND_MOVES(Type) \
  NOTHROW_DEFAULT_CTOR(Type);           \
  NOTHROW_MOVES(Type)

namespace pplite {

using dim_type = int;
using size_t = std::size_t;

enum class Topol { CLOSED, NNC };
enum class Spec_Elem { EMPTY, UNIVERSE };
enum class Widen_Spec { SAFE, RISKY };
enum class Widen_Impl { H79, BOXED_H79, BHRZ03 };

constexpr dim_type not_a_dim() { return -1; }

template <typename T>
constexpr dim_type sizeof_to_bits(const T size) {
  return size * std::numeric_limits<unsigned char>::digits;
}

using Dims = std::vector<dim_type>;

enum class TV_Bool { DONT_KNOW, FALSE, TRUE };

using Counter = unsigned long long;

namespace detail{

// If not set, default topology is CLOSED.
extern PPLITE_TLS Topol default_topol;

// If not set, widening specification is SAFE.
extern PPLITE_TLS Widen_Spec widen_spec;

// If not set, widening implementation is H79.
extern PPLITE_TLS Widen_Impl widen_impl;

} // namespace detail

inline Topol get_default_topology() {
  return detail::default_topol;
}
inline void set_default_topology(Topol topol) {
  detail::default_topol = topol;
}

inline Widen_Spec get_widen_spec() {
  return detail::widen_spec;
}
inline void set_widen_spec(Widen_Spec w_spec) {
  detail::widen_spec = w_spec;
}

inline Widen_Impl get_widen_impl() {
  return detail::widen_impl;
}
inline void set_widen_impl(Widen_Impl w_impl) {
  detail::widen_impl = w_impl;
}

} // namespace pplite

#endif // !defined(pplite_globals_hh)
