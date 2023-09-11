/* PPLite: a lightweight library for convex polyhedra derived from PPL.
   Copyright (C) 2001-2010 Roberto <bagnara@cs.unipr.it>
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

#include "pplite_test.hh"
#include <vector>
#include <functional>
#include <cmath>

#ifndef M_PI
# define M_PI           3.14159265358979323846  /* pi */
#endif

int
count_points(const Poly& ph) {
  assert(ph.is_necessarily_closed());
  ph.minimize();
  const auto& gs = ph.gens();
  return std::count_if(gs.begin(), gs.end(), std::mem_fn(&Gen::is_point));
}

bool
test01() {
  // Set up a random numbers' generator.
  gmp_randclass rg(gmp_randinit_default);

  Var x(0);
  Var y(1);
  Var z(2);

  const int maxc = 10000;
  Poly ph(3, Spec_Elem::EMPTY);
  nout << count_points(ph) << endl;
  for (int n = 0; n < 200; ++n) {
    const Integer cx = mpz_class(rg.get_z_range(maxc)).get_mpz_t();
    const Integer cy = mpz_class(rg.get_z_range(maxc)).get_mpz_t();
    const Integer cz = mpz_class(rg.get_z_range(maxc)).get_mpz_t();
    ph.add_gen(point(cx*x + cy*y + cz*z));
    if (ph.is_empty())
      return false;
    nout << count_points(ph) << endl;
  }
  return true;
}

void
point_on_the_unit_n_sphere(unsigned n,
                           const std::vector<float>& theta,
                           std::vector<float>& coordinate) {
  assert(n >= 2);

  if (n == 2) {
    coordinate[0] *= sin(theta[0]);
    coordinate[1] *= cos(theta[0]);
  }
  else {
    point_on_the_unit_n_sphere(n-1, theta, coordinate);
    float sin_theta_n_2 = sin(theta[n-2]);
    for (unsigned i = n-1; i-- > 0; )
      coordinate[i] *= sin_theta_n_2;
    coordinate[n-1] *= cos(theta[n-2]);
  }
}

void
random_polytope(Poly& ph, unsigned dim, unsigned num_points,
                float radius = 1.0) {
  assert(dim >= 2);

  std::vector<float> theta(dim - 1);
  std::vector<float> coordinate(dim);

  for (unsigned n = num_points; n > 0; --n) {
    // Compute n-1 random angles.
    for (unsigned i = dim - 1; i-- > 0; )
      theta[i] = 2.0*M_PI*static_cast<double>(rand())/RAND_MAX;
    // Compute float coordinates.
    for (unsigned i = dim; i-- > 0; )
      coordinate[i] = radius;
    point_on_the_unit_n_sphere(dim, theta, coordinate);

    Linear_Expr le;
    for (unsigned i = dim; i-- > 0; ) {
      mpz_class z = coordinate[i]*1000000.0;
      add_mul_assign(le, Integer(z.get_mpz_t()), Var(i));
    }
    ph.add_gen(point(le));
  }
}

bool
test02() {
#ifdef NDEBUG
  const dim_type limit = 7;
#else
  const dim_type limit = 4;
#endif

  for (dim_type dim = 2; dim <= limit; ++dim) {
    Poly ph(dim, Spec_Elem::EMPTY);
    random_polytope(ph, dim, dim * dim);
    ph.minimize();
    Gens gs = ph.copy_gens();
    if (!std::all_of(gs.begin(), gs.end(), std::mem_fn(&Gen::is_point)))
      return false;

    nout << "dimension = " << dim
         << ", points = " << ph.num_min_gens() << " (" << (dim * dim) << ")"
         << ", constraints = " << ph.num_min_cons() << endl;
  }
  return true;
}

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
END_MAIN
