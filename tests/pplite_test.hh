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

#ifndef PPLite_pplite_test_hh
#define PPLite_pplite_test_hh 1

#include "pplite.hh"

#include <stdexcept>
#include <sstream>
#include <list>
#include <iterator>
#include <string>
#include <iostream>
#include <algorithm>
#include <typeinfo>
#include <cstdlib>

// For debugging purposes (only meaningful with assertions on).
/* #define Poly Two_Poly<B_Poly, Poly> */
/* #define Poly Two_Poly<F_Poly, Poly> */

#define MyPSet P_Set
#define MyDisj MyPSet::Disj

// Macro allowing for reusing Poly tests on Dyn_Poly.
// Example usage: compile with flag -DTEST_ABS_POLY_KIND=U_POLY
#ifdef TEST_ABS_POLY_KIND
#define Poly dynamic::Dyn_Poly
#else
#define TEST_ABS_POLY_KIND POLY
#endif

#define BEGIN_MAIN                                        \
int main() try {                                          \
  using namespace dynamic;                                \
  auto kind = Abs_Poly::Kind::TEST_ABS_POLY_KIND;         \
  auto kind_name = abs_poly_kind_to_name(kind);           \
  if (not set_default_poly_kind(kind_name, false)) {      \
    std::cerr << "Failed setting poly kind: ";            \
    std::cerr << kind_name << "\n";                       \
    return 1;                                             \
  }                                                       \
  bool succeeded = false;                                 \
  std::list<std::string> failed_tests;                    \
  std::list<std::string> unexpectedly_succeeded_tests;

#define END_MAIN                                                        \
  if (!failed_tests.empty()) {                                          \
    std::cerr << "tests failed: ";                                      \
    std::copy(failed_tests.begin(),                                     \
              failed_tests.end(),                                       \
              std::ostream_iterator<std::string>(std::cerr, " "));      \
    std::cerr << std::endl;                                             \
    return 1;                                                           \
  }                                                                     \
  if (!unexpectedly_succeeded_tests.empty()) {                          \
    std::cerr << "tests unexpectedly succeeded: ";                      \
    std::copy(unexpectedly_succeeded_tests.begin(),                     \
              unexpectedly_succeeded_tests.end(),                       \
              std::ostream_iterator<std::string>(std::cerr, " "));      \
    std::cerr << std::endl;                                             \
    return 1;                                                           \
  }                                                                     \
  return 0;                                                             \
}                                                                       \
catch (const std::exception& e) {                                       \
  std::cerr << "std::exception caught: "                                \
            << e.what() << " (type == " << typeid(e).name() << ")"      \
            << std::endl;                                               \
  exit(1);                                                              \
}

#define ANNOUNCE_TEST(test)              \
  nout << "\n=== " #test " ===" << std::endl

#define RUN_TEST(test)                                                  \
  try {                                                                 \
    succeeded = test();                                                 \
  }                                                                     \
  catch (const std::exception& e) {                                     \
    nout << "std::exception caught: "                                   \
         << e.what() << " (type == " << typeid(e).name() << ")"         \
         << std::endl;                                                  \
    succeeded = false;                                                  \
  }                                                                     \
  catch (...) {                                                         \
    nout << "unknown exception caught" << std::endl;                    \
    succeeded = false;                                                  \
  }

#define DO_TEST(test)                    \
  ANNOUNCE_TEST(test);                   \
  RUN_TEST(test);                        \
  if (!succeeded)                        \
    failed_tests.push_back(#test);

#define DO_TEST_F(test)                                 \
  ANNOUNCE_TEST(test);                                  \
  RUN_TEST(test);                                       \
  if (succeeded)                                        \
    unexpectedly_succeeded_tests.push_back(#test);

namespace pplite {
namespace test {

inline bool
check_exp_eval() {
  return getenv("PPLITE_EXP_EVAL") != 0;
}

inline bool
check_noisy(const char* environment_variable) {
  return getenv(environment_variable) != 0;
}

template<typename CharT, typename Traits = std::char_traits<CharT> >
class nullbuf : public std::basic_streambuf<CharT, Traits> {
protected:
  virtual typename Traits::int_type overflow(typename Traits::int_type c) {
    return Traits::not_eof(c);
  }
};

template <class CharT, class Traits = std::char_traits<CharT> >
class noisy_ostream : public std::basic_ostream<CharT, Traits> {
private:
  nullbuf<CharT, Traits> black_hole;

public:
  noisy_ostream(const std::basic_ostream<CharT, Traits>& os,
                const char* environment_variable)
    : std::basic_ostream<CharT, Traits>(check_noisy(environment_variable)
                                        ? os.rdbuf()
                                        : &black_hole) {
  }
};

static noisy_ostream<char> nout(std::cout, "PPLITE_NOISY_TESTS");
static noisy_ostream<char> vnout(std::cout, "PPLITE_VERY_NOISY_TESTS");

inline bool
check_print(const std::string& s1, const std::string& s2) {
  if (s1 == s2)
    return true;
  nout << "\n=== Printing mismatch ===\n";
  nout << "Printing lhs:\n";
  nout << s1;
  nout << "\nPrinting rhs:\n";
  nout << s2;
  nout << "\n=========================\n";
  return false;
}

template <typename T>
bool check_print(const T& t, const std::string& s) {
  std::stringstream sst;
  t.print(sst);
  return check_print(sst.str(), s);
}

template <typename T, typename U>
bool print_same(const T& t, const U& u) {
  std::stringstream sst;
  t.print(sst);
  std::stringstream ssu;
  u.print(ssu);
  return check_print(sst.str(), ssu.str());
}

template <typename Iter>
void print_seq(Iter first, Iter last, const char* sep = "\n") {
  using IO_Operators::operator<<;
  for ( ; first != last; ++first)
    nout << *first << sep;
}

template <typename T>
void print(const std::vector<T>& ts) { print_seq(ts.cbegin(), ts.cend()); }

inline void
print_cons(const Cons& cs) { print(cs); }
inline void
print_gens(const Gens& gs) { print(gs); }

template <typename PH>
inline void
print_cons(const PH& ph, const std::string& s) {
  nout << s << std::endl;
  const auto& cs = ph.cons();
  print_seq(cs.cbegin(), cs.cend());
}

template <typename PH>
inline void
print_gens(const PH& ph, const std::string& s) {
  nout << s << std::endl;
  const auto& gs = ph.gens();
  print_seq(gs.cbegin(), gs.cend());
}

template <typename PH>
inline void
print_pset(const PolySet<PH>& ps, const std::string& s) {
  nout << s << std::endl;
  dim_type i = 0;
  for (const auto& ph : ps) {
    nout << "== disj " << i << " ===\n";
    ++i;
    const auto& cs = ph.cons();
    print_seq(cs.cbegin(), cs.cend());
  }
}

} // namespace test
} // namespace pplite

// These using directive and declaration are just to avoid the
// corresponding namespace qualifications in all the tests.
using namespace pplite;
using namespace pplite::test;
using std::endl;

#endif // !defined(PPLite_pplite_test_hh)
