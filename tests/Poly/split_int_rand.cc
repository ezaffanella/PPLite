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

#include "pplite_test.hh"

#include <random>

using namespace IO_Operators;

// Avoid seeding (for repeatability).
// std::random_device ran_dev;
// std::mt19937 rng(ran_dev());
std::mt19937 rng;
const int abs_max = 100000;
std::uniform_int_distribution<int> myrand(-abs_max, abs_max);

Con
random_homo_con(dim_type dim) {
  assert(dim >= 2);
  Linear_Expr le;
  for (dim_type i = dim; i-- > 0; ) {
    Integer coeff = myrand(rng);
    add_mul_assign(le, coeff, Var(i));
  }
  return le >= 0;
}

Cons
perturbed_hypercube(dim_type dim) {
  Linear_Expr pert;
  for (dim_type i = 0; i < dim; ++i) {
    dim_type z = ((i % 2 == 0) ? -1 : 1) * (i+7);
    add_mul_assign(pert, Integer(z), Var(i));
  }

  Cons cs;
  cs.reserve(2 * dim);
  for (dim_type i = 0; i < dim; ++i) {
    Linear_Expr le = pert;
    add_mul_assign(le, Integer(abs_max), Var(i));
    cs.push_back(le <= abs_max);
    cs.push_back(le >= -abs_max);
  }
  return cs;
}

Poly get_hypercube(dim_type dim) {
  Poly ph(dim, Topol::CLOSED);
  ph.add_cons(perturbed_hypercube(dim));
  ph.minimize();
  nout << "\n=== dimension = " << dim
       << ", gens = " << ph.num_min_gens()
       << ", cons = " << ph.num_min_cons()
       << " ===" << endl;
  return ph;
}

void
new_split(const Poly& ph, const Con& c1) {
  nout << "\n=== Time for new style integral split: ";
  Clock clock;
  auto ph1 = ph;
  auto ph2 = ph1.integral_split(c1);
  clock.print_elapsed(nout);
  nout << endl;
  nout << "pos_split: "
       << "gens = " << ph1.num_min_gens() << ", "
       << "cons = " << ph1.num_min_cons() << endl;
  nout << "neg_split: "
       << "gens = " << ph2.num_min_gens() << ", "
       << "cons = " << ph2.num_min_cons() << endl;
}

void old_split(const Poly& ph, const Con& c) {
  nout << "\n=== Time for old style integral split: ";
  Clock clock;
  auto ph1 = ph;
  auto [c1, c2] = detail::integral_complement_cons(c);
  auto ph2 = ph1;
  ph1.add_con(c1);
  ph2.add_con(c2);
  ph1.minimize();
  ph2.minimize();
  clock.print_elapsed(nout);
  nout << endl;
  nout << "pos_split: "
       << "gens = " << ph1.num_min_gens() << ", "
       << "cons = " << ph1.num_min_cons() << endl;
  nout << "neg_split: "
       << "gens = " << ph2.num_min_gens() << ", "
       << "cons = " << ph2.num_min_cons() << endl;
}

bool
test01() {
  dim_type min_dim = 5;
  dim_type max_dim = 7;
  if (check_exp_eval()) {
    min_dim = 10;
    max_dim = 12;
  }

  for (dim_type dim = min_dim; dim <= max_dim; ++dim) {
    Poly ph = get_hypercube(dim);
    Con con = random_homo_con(dim);
    nout << "Splitting on constraint: " << con << endl;
    old_split(ph, con);
    new_split(ph, con);
  }

  return true;
}

bool
test02() {
  dim_type min_dim = 5;
  dim_type max_dim = 7;
  if (check_exp_eval()) {
    min_dim = 12;
    max_dim = 14;
  }

  using namespace IO_Operators;

  for (dim_type dim = min_dim; dim <= max_dim; ++dim) {
    Poly ph(dim, Topol::CLOSED);
    ph.add_cons(perturbed_hypercube(dim));
    ph.minimize();
    nout << "\n=== dimension = " << dim
         << ", gens = " << ph.num_min_gens() << " (" << (dim * dim) << ")"
         << ", cons = " << ph.num_min_cons() << " ===" << endl;

    Con c = random_homo_con(dim);
    nout << "Splitting on constraint: " << c << endl;

    Poly ph1;
    Poly ph2;
    {
      nout << "\n=== Time for old style integral split: ";
      Clock clock;
      ph1 = ph;
      ph2 = ph;
      auto [c1, c2] = detail::integral_complement_cons(c);
      ph1.add_con(c1);
      ph2.add_con(c2);
      ph1.minimize();
      ph2.minimize();
      clock.print_elapsed(nout);
    }
    {
      nout << "\n=== Time for new style integral split: ";
      Clock clock;
      auto ph1_split = ph;
      auto ph2_split = ph1_split.integral_split(c);
      clock.print_elapsed(nout);

      // Check:
      bool ok = (ph1 == ph1_split && ph2 == ph2_split);
      if (!ok) {
        nout << "Difference detected!\n";
        return false;
      }
    }
    nout << endl;
  }
  return true;
}

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
END_MAIN
