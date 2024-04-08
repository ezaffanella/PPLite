/* PPLite: a lightweight library for convex polyhedra derived from PPL.
   Copyright (C) 2001-2010 Roberto <bagnara@cs.unipr.it>
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

#include "pplite_test.hh"
#include <cmath>

namespace {

Con
complement_con(const Con& beta, bool nnc=true) {
  auto expr = beta.linear_expr();
  auto inhomo = beta.inhomo_term();
  neg_assign(inhomo);
  neg_assign(expr);
  Con::Type type = Con::NONSTRICT_INEQUALITY;
  if (nnc && beta.is_nonstrict_inequality())
      type = Con::STRICT_INEQUALITY;
  return Con(expr, inhomo, type);
}

#define LAZY_SPLIT 0

bool
test01() {
  const dim_type first_dim = 4;
  const dim_type last_dim = 4;
  const dim_type split_nr = 4;
  bool ok = true;

  for (dim_type dim = first_dim; dim <= last_dim; ++dim) {
    const dim_type split_dim = dim;
    Cons hyper_cs;
    hyper_cs.reserve(2 * dim);
    for (dim_type i = 0; i < dim; ++i) {
      if (i%2) {
        hyper_cs.push_back(Var(i) > 0);
        hyper_cs.push_back(Var(i) < 10);
      } else {
        hyper_cs.push_back(Var(i) >= 0);
        hyper_cs.push_back(Var(i) <= 10);
      }
    }

    Poly ph(dim, Topol::NNC);
    Cons cons;
    cons.reserve(split_nr*split_dim);
    for (dim_type i = 0; i < split_dim; ++i)
      for (dim_type j = 1; j <= split_nr; ++j)
        cons.push_back(Var(i) < j);

    std::vector<Poly> ph_split;
    std::vector<Poly> knres_eager;
    const dim_type max_nr = std::pow((split_nr + 1), (split_dim + 1));
    ph_split.resize(max_nr, ph);
    knres_eager.resize(max_nr, ph);
#if LAZY_SPLIT
    std::vector<Poly> knres_lazy;
    knres_lazy.resize(max_nr, ph);
#endif

    // Build hypercube
    ph.add_cons(hyper_cs);
    ph.minimize();
    knres_eager[0] = ph;
    ph_split[0] = ph;
#if LAZY_SPLIT
    knres_lazy[0] = ph;
#endif

    nout << "\n*************************\n"
         << "Hypercube of dim " << dim
         << " split " << split_nr << " times on each dim.\n";

    nout << "\nNew split mode: \t | ";
    {
      dim_type c = 0;
      Clock clock;
      for (dim_type k = 0; k < split_dim; ++k) {
        const dim_type step = std::pow((split_nr + 1), k);
        for (dim_type n = 0; n < split_nr; ++n) {
          for (dim_type i = 0; i < step; ++i) {
            ph_split[i + (n+1)*step] =
              ph_split[i + n * step].split(cons[c]);
          }
          ++c;
        }
      }
      clock.print_elapsed(nout);
    }

    nout << "\nUsual eager mode \t | ";
    {
      dim_type c = 0;
      Clock clock;
      for (dim_type k = 0; k < split_dim; ++k) {
        const dim_type step = std::pow((split_nr + 1), k);
        for (dim_type n = 0; n < split_nr; ++n) {
          for (dim_type i = 0; i < step; ++i) {
            knres_eager[i + (n+1)*step] = knres_eager[i + n * step];
            knres_eager[i + n * step].add_con(cons[c]);
            knres_eager[i + n * step].minimize();
            knres_eager[i + (n+1)*step].add_con(complement_con(cons[c]));
            knres_eager[i + (n+1)*step].minimize();
          }
          ++c;
        }
      }
      clock.print_elapsed(nout);
    }

    // Check eager
    for (dim_type i = 0; i < max_nr; ++i) {
      ok = (ph_split[i] == knres_eager[i]);
      if (!ok) {
        nout <<"\n Difference detected! \n";
        nout << "****** ph_then ********* ";
        ph_split[i].ascii_dump(nout);
        nout << "****** knres_then ********* ";
        knres_eager[i].ascii_dump(nout);
        return ok;
      }
    }

#if LAZY_SPLIT
    nout << "\nUsual lazy mode \t | ";
    {
      dim_type c = 0;
      Clock clock;
      for (dim_type k = 0; k < split_dim; ++k) {
        const dim_type step = std::pow((split_nr + 1), k);
        for (dim_type n = 0; n < split_nr; ++n) {
          for (dim_type i = 0; i < step; ++i) {
            knres_lazy[i + (n+1)*step] = knres_lazy[i + n * step];
            knres_lazy[i + n * step].add_con(cons[c]);
            knres_lazy[i + (n+1)*step].add_con(complement_con(cons[c]));
          }
          ++c;
        }
      }
      for (dim_type i = 0; i < max_nr; ++i)
        knres_lazy[i].minimize();

      clock.print_elapsed(nout);
    }
    // Check lazy
    for (dim_type i = 0; i < max_nr; ++i) {
      ok = (ph_split[i] == knres_lazy[i]);
      if (!ok) {
        nout <<"\n Difference detected! \n";
        nout << "****** ph_then ********* ";
        ph_split[i].ascii_dump(nout);
        nout << "****** knres_then ********* ";
        knres_lazy[i].ascii_dump(nout);
        return ok;
      }
    }
#endif
    nout << "\n";
  }
  return ok;
}

} // namespace

BEGIN_MAIN
  DO_TEST(test01);
END_MAIN
