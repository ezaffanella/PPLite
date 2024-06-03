/* PPLite: a lightweight library for convex polyhedra derived from PPL.
   Copyright (C) 2001-2010 Roberto Bagnara <bagnara@cs.unipr.it>
   Copyright (C) 2010-2018 BUGSENG srl (http://bugseng.com)
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

#include "pplite.hh"

#include <string>
#include <vector>
#include <set>
#include <limits>
#include <climits>
#include <cassert>
#include <cstdarg>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <stdexcept>

#include <getopt.h>

using namespace pplite;

using poly_type = pplite::Poly;
// using poly_type = pplite::Poly_Stats;

namespace {

const char* program_name = nullptr;
bool print_timings = false;
bool sort_input = false;
bool verbose = false;
const char* check_file_name = nullptr;

const char* input_file_name = nullptr;
std::istream* input_stream_p = nullptr;
const char* output_file_name = nullptr;
std::ostream* output_stream_p = nullptr;

const char* const option_letters = "ho:stvc:";

struct option long_options[] = {
  {"help",    no_argument,       nullptr, 'h'},
  {"output",  required_argument, nullptr, 'o'},
  {"sort-input", no_argument,    nullptr, 's'},
  {"timings", no_argument,       nullptr, 't'},
  {"verbose", no_argument,       nullptr, 'v'},
  {"check",   required_argument, nullptr, 'c'},
  {nullptr, 0, nullptr, 0}
};

void
usage(const std::string& program) {
  std::cout
    << "Usage: " << program << " [OPTION]... [FILE]\n"
    "Reads an H-representation (resp., a V-representation) of a polyhedron\n"
    "and generates a V-representation (resp., an H-representation) of\n"
    "the same polyhedron.\n\n"
    "Options:\n"
    "  -h, --help              prints this help text to stdout\n"
    "  -oPATH, --output=PATH   appends output to PATH\n"
    "  -t, --timings           prints timings to stderr\n"
    "  -v, --verbose           produces lots of output\n"
    "  -cPATH, --check=PATH    checks if result is equal to what is in PATH\n";
}

[[ noreturn ]] void
fatal(const std::string& msg) {
  std::cerr << "Fatal error: " << msg << "\n";
  exit(1);
}

[[ noreturn ]] void
error(const std::string& msg) {
  std::cerr << "Error: " << msg << "\n";
  exit(1);
}

void
warning(const std::string& msg) {
  std::cerr << "Warning: " << msg << "\n";
}

void
set_input(const char* const file_name) {
  if (input_stream_p != &std::cin)
    delete input_stream_p;

  if (file_name) {
    input_stream_p = new std::ifstream(file_name, std::ios_base::in);
    if (!*input_stream_p)
      fatal(std::string("cannot open input file `") + file_name + "'");
    input_file_name = file_name;
  } else {
    input_stream_p = &std::cin;
    input_file_name = "<cin>";
  }
}

std::istream&
input() {
  assert(input_stream_p != nullptr);
  return *input_stream_p;
}

void
set_output(const char* const file_name) {
  if (output_stream_p != &std::cout)
    delete output_stream_p;

  if (file_name) {
    output_stream_p = new std::ofstream(file_name,
                                        std::ios_base::out
                                        | std::ios_base::app);
    if (!*output_stream_p)
      fatal(std::string("cannot open output file `") + file_name + "'");
    output_file_name = file_name;
  } else {
    output_stream_p = &std::cout;
    output_file_name = "<cout>";
  }
}

std::ostream&
output() {
  assert(output_stream_p != 0);
  return *output_stream_p;
}

void
process_options(const int argc, char* const argv[]) {
  while (true) {
    int option_index = 0;
    const int c = getopt_long(argc, argv,
                              option_letters, long_options,
                              &option_index);
    if (c == EOF)
      break;

    switch (c) {
    case 0:
      break;

    case '?':
    case 'h':
      usage(argv[0]);
      exit(0);
      break;

    case 'o':
      output_file_name = optarg;
      break;

    case 's':
      sort_input = true;
      break;

    case 't':
      print_timings = true;
      break;

    case 'v':
      verbose = true;
      break;

    case 'c':
      check_file_name = optarg;
      break;

    default:
      abort();
    }
  }

  const auto num_input_files = argc - optind;
  if (num_input_files > 1)
    fatal("at most one input file is accepted");
  if (num_input_files == 1)
    input_file_name = argv[optind];
}

void
maybe_print_clock(const Clock& clock) {
  if (print_timings) {
    std::cerr << input_file_name << " ";
    clock.print_elapsed(std::cerr);
    std::cerr << std::endl;
  }
}

template <typename T>
bool
guarded_read(std::istream& in, T& x) {
  try {
    in >> x;
    return !in.fail();
  }
  catch (...) {
    return false;
  }
}

template <typename T>
inline void
guarded_read(std::istream& in, T& x, const std::string& error_msg) {
  if (!guarded_read(in, x))
    error(error_msg);
}

template <typename T>
void
guarded_write(std::ostream& out, const T& x) {
  bool succeeded = false;
  try {
    out << x;
    succeeded = !out.fail();
  }
  catch (...) {
  }
  if (!succeeded) {
    fatal(std::string("cannot write to output file `")
          + output_file_name + "'");
  }
}

enum class Num_Type { INTEGER, RATIONAL, REAL };

void
read_coefficients(std::istream& in,
                  const Num_Type num_type,
                  Integers& coeffs,
                  Integer& denom) {
  const auto num_coeffs = coeffs.size();
  switch (num_type) {
  case Num_Type::INTEGER:
    {
      for (unsigned i = 0; i < num_coeffs; ++i)
        guarded_read(in, coeffs[i],
                     "missing or invalid integer coefficient");
      denom = 1;
      break;
    }
  case Num_Type::RATIONAL:
    {
      Rationals rational_coeffs(num_coeffs);
      for (unsigned i = 0; i < num_coeffs; ++i)
        guarded_read(in, rational_coeffs[i],
                     "missing or invalid rational coefficient");
      rationals_to_integers(rational_coeffs, coeffs, denom);
      break;
    }
  case Num_Type::REAL:
    {
      Rationals rational_coeffs;
      rational_coeffs.reserve(num_coeffs);
      for (unsigned i = 0; i < num_coeffs; ++i) {
        double d;
        guarded_read(in, d, "missing or invalid real coefficient");
        rational_coeffs.emplace_back(d);
      }
      rationals_to_integers(rational_coeffs, coeffs, denom);
      break;
    }
  }
}

void
read_indexes_set(std::istream& in, std::set<unsigned>& dst) {
  assert(dst.empty());
  unsigned num_elems;
  guarded_read(in, num_elems, "missing or invalid number of set elements");
  while (num_elems--) {
    unsigned i;
    guarded_read(in, i, "missing or invalid set element");
    dst.insert(i);
  }
}

template <typename Row>
struct RowCmp {
  bool operator()(const Row& row1, const Row& row2) const {
    return compare(row1, row2) < 0;
  }
};

enum class Repr { H, V };

Repr
read_polyhedron(std::istream& in, poly_type& ph) {
  // By default we have an H-representation.
  Repr repr = Repr::H;

  std::string s;
  std::set<unsigned> linearity;
  while (true) {
    guarded_read(in, s, "premature end of file while seeking for `begin'");
    if (s == "V-representation")
      repr = Repr::V;
    else if (s == "H-representation")
      repr = Repr::H;
    else if (s == "linearity" || s == "equality" || s == "partial_enum") {
      read_indexes_set(in, linearity);
      if (verbose) {
        std::cerr << "Linearity: ";
        for (auto i : linearity)
          std::cerr << i << " ";
        std::cerr << std::endl;
      }
    }
    else if (s == "begin")
      break;
    else
      // A comment: skip to end of line.
      in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }

  unsigned num_rows;
  guarded_read(in, num_rows, "illegal or missing number of rows");
  unsigned num_cols;
  guarded_read(in, num_cols, "illegal or missing number of columns");
  assert(num_cols > 0);
  const dim_type space_dim = num_cols - 1;

  guarded_read(in, s, "missing number type");
  Num_Type num_type = Num_Type::INTEGER;
  if (s == "integer")
    num_type = Num_Type::INTEGER;
  else if (s == "rational")
    num_type = Num_Type::RATIONAL;
  else if (s == "real")
    num_type = Num_Type::REAL;
  else
    error(std::string("illegal number type `" + s + "'"));

  if (verbose) {
    std::cerr << "Problem dimension: "
              << num_rows << " x " << num_cols
              << "; number type: " << s
              << std::endl;
  }

  Cons cs;
  Gens gs;

  unsigned row = 0;
  const auto linearity_end = linearity.end();
  if (repr == Repr::V) {
    // The V representation allows for `space_dim' coordinates.
    Integers coeffs(space_dim);
    Integer denom;
    bool has_a_point = false;
    for (row = 0; row < num_rows; ++row) {
      int vertex_marker;
      guarded_read(in, vertex_marker, "missing vertex marker");
      if (vertex_marker < 0 || vertex_marker > 1)
        error("illegal vertex marker");
      read_coefficients(in, num_type, coeffs, denom);
      Linear_Expr e;
      for (dim_type j = space_dim; j-- > 0; )
        add_mul_assign(e, coeffs[j], Var(j));
      if (vertex_marker == 1) {
        assert(linearity.find(row+1) == linearity_end);
        gs.push_back(point(std::move(e), std::move(denom)));
        has_a_point = true;
      }
      else if (linearity.find(row+1) != linearity_end)
        gs.push_back(line(std::move(e)));
      else
        gs.push_back(ray(std::move(e)));
    }
    // Every non-empty generator system must have at least one point.
    if (row > 0 && !has_a_point)
      gs.push_back(point());
    if (verbose) {
      using namespace pplite::IO_Operators;
      std::cerr << "Generator system:\n" << gs << std::endl;
    }
  }
  else {
    assert(repr == Repr::H);
    // The H representation stores the inhomogeneous term at index 0,
    // and the variables' coefficients at indices 1, 2, ..., space_dim.
    Integers coeffs(space_dim+1);
    Integer denom;
    for (row = 0; row < num_rows; ++row) {
      read_coefficients(in, num_type, coeffs, denom);
      Linear_Expr e;
      for (dim_type j = num_cols; j-- > 1; )
        add_mul_assign(e, coeffs[j], Var(j-1));
      neg_assign(coeffs[0]);
      if (linearity.find(row+1) != linearity_end)
        cs.push_back(e == coeffs[0]);
      else
        cs.push_back(e >= coeffs[0]);
    }
    if (verbose) {
      using namespace pplite::IO_Operators;
      std::cerr << "Constraint system:\n" << cs << std::endl;
    }
  }

  guarded_read(in, s, "premature end of file while seeking for `end'");
  if (s != "end")
    error(std::string("found `") + s + " while seeking for `end'");

  if (repr == Repr::H) {
    if (sort_input)
      std::sort(cs.begin(), cs.end(), RowCmp<Con>());
    ph = poly_type(space_dim);
    ph.add_cons(std::move(cs));
  } else {
    if (sort_input)
      std::sort(gs.begin(), gs.end(), RowCmp<Gen>());
    ph = poly_type(space_dim, Spec_Elem::EMPTY);
    ph.add_gens(std::move(gs));
  }
  return repr;
}

void
write_polyhedron(std::ostream& out, const poly_type& ph, const Repr repr) {
  Cons cons;
  Gens gens;
  dim_type num_rows = 0;
  if (repr == Repr::H) {
    cons = ph.copy_cons();
    num_rows = cons.size();
    if (verbose) {
      using namespace pplite::IO_Operators;
      std::cerr << "Constraint system:\n" << cons << std::endl;
    }
  } else {
    assert(repr == Repr::V);
    gens = ph.copy_gens();
    num_rows = gens.size();
    if (verbose) {
      using namespace pplite::IO_Operators;
      std::cerr << "Generator system:\n" << gens << std::endl;
    }
  }

  guarded_write(out, (repr == Repr::H) ? "H" : "V");
  guarded_write(out, "-representation\n");

  std::vector<unsigned> linearity;
  if (repr == Repr::H) {
    for (auto i = 0; i < num_rows; ++i) {
      if (cons[i].is_equality())
        linearity.push_back(i);
    }
  } else {
    for (auto i = 0; i < num_rows; ++i) {
      if (gens[i].is_line())
        linearity.push_back(i);
    }
  }

  if (!linearity.empty()) {
    guarded_write(out, "linearity ");
    guarded_write(out, linearity.size());
    for (auto i : linearity) {
      guarded_write(out, ' ');
      // Note: 1-based indexing.
      guarded_write(out, i+1);
    }
    guarded_write(out, '\n');
  }

  const auto space_dim = ph.space_dim();

  guarded_write(out, "begin\n");
  guarded_write(out, num_rows);
  guarded_write(out, ' ');
  guarded_write(out, space_dim + 1);
  guarded_write(out, ' ');
  if (repr == Repr::H)
    guarded_write(out, "integer\n");
  else
    guarded_write(out, "rational\n");

  if (repr == Repr::H) {
    for (const auto& c : cons) {
      guarded_write(out, c.inhomo_term());
      for (dim_type j = 0; j < space_dim; ++j) {
        guarded_write(out, ' ');
        guarded_write(out, c.coeff(Var(j)));
      }
      guarded_write(out, '\n');
    }
  }
  else {
    assert(repr == Repr::V);
    for (const auto& g : gens) {
      if (g.is_point()) {
        guarded_write(out, '1');
        const auto& divisor = g.divisor();
        for (dim_type j = 0; j < space_dim; ++j) {
          guarded_write(out, ' ');
          if (g.coeff(Var(j)) == 0)
            guarded_write(out, '0');
          else {
            const auto& numer = g.coeff(Var(j));
            guarded_write(out, Rational(numer, divisor));
          }
        }
      } else {
        assert(g.is_ray() || g.is_line());
        guarded_write(out, '0');
        for (dim_type j = 0; j < space_dim; ++j) {
          guarded_write(out, ' ');
          guarded_write(out, g.coeff(Var(j)));
        }
      }
      guarded_write(out, '\n');
    }
  }
  guarded_write(out, "end\n");

  // Flush `out'.
  bool flush_succeeded = false;
  try {
    out.flush();
    flush_succeeded = !out.fail();
  }
  catch (...) {
  }
  if (!flush_succeeded) {
    fatal(std::string("cannot write to output file `")
          + output_file_name + "'");
  }
}

#define CATCH_ALL               \
catch (const std::bad_alloc&) { \
  fatal("out of memory");       \
  exit(1);                      \
}                               \
catch (...) {                   \
  fatal("internal error");      \
  exit(1);                      \
}

void
convert(const char* name_in) try {
  // Set up the input and output streams.
  set_input(name_in);
  set_output(output_file_name);

  poly_type ph;
  const Repr repr = read_polyhedron(input(), ph);

  // Warn for misplaced linearity commands, and ignore all what follows.
  std::string s;
  while (guarded_read(input(), s)) {
    if (s == "linearity" || s == "equality" || s == "partial_enum") {
      error("the `linearity' command must occur before `begin'");
    }
    input().ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }

  Clock clock;

  enum Command { None, H_to_V, V_to_H };
  Command command = (repr == Repr::V) ? V_to_H : H_to_V;

  // Compute the dual representation.
  ph.minimize();

  maybe_print_clock(clock);

  // Write the result of the conversion.
  if (repr == Repr::V)
    write_polyhedron(output(), ph, Repr::H);
  else
    write_polyhedron(output(), ph, Repr::V);

  // Check the result, if requested to do so.
  if (check_file_name) {
    set_input(check_file_name);
    // Read the polyhedron containing the expected result.
    poly_type e_ph;
    const Repr e_repr = read_polyhedron(input(), e_ph);

    switch (command) {
    case H_to_V:
      {
        if (e_repr == Repr::H)
          warning("checking an H-to-V conversion with an H representation");
        // If the polyhedra differ, that is the problem.
        if (ph != e_ph) {
          if (verbose)
            std::cerr << "Check failed: polyhedra differ" << std::endl;
          exit(1);
        }
        const auto ph_num_gens = ph.num_min_gens();
        const auto e_ph_num_gens = e_ph.num_min_gens();
        if (ph_num_gens != e_ph_num_gens) {
          // If we have different number of generators, we fail.
          std::cerr << "Check failed: different number of generators:\n"
                    << "expected " << e_ph_num_gens
                    << ", obtained " << ph_num_gens
                    << std::endl;
          exit(1);
        }
        break;
      }
    case V_to_H:
      {
        if (e_repr == Repr::V)
          warning("checking an V-to-H conversion with a V representation");
        // If the polyhedra differ, that is the problem.
        if (ph != e_ph) {
          if (verbose)
            std::cerr << "Check failed: polyhedra differ" << std::endl;
          exit(1);
        }
        const auto ph_num_cons = ph.num_min_cons();
        const auto e_ph_num_cons = e_ph.num_min_cons();
        if (ph_num_cons != e_ph_num_cons) {
          // If we have different number of constraints, we fail.
          std::cerr << "Check failed: different number of constraints:\n"
                    << "expected " << e_ph_num_cons
                    << ", obtained " << ph_num_cons
                    << std::endl;
          exit(1);
        }
        break;
      }
    case None:
      break;
    }
  }

}
CATCH_ALL

} // namespace

int
main(int argc, char* argv[]) try {
  program_name = argv[0];
  // Process command line options.
  process_options(argc, argv);
  // Passing nullptr means "read from stdin".
  convert(input_file_name);
  return 0;
}
CATCH_ALL
