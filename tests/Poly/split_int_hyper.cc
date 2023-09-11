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

namespace {

template <typename Iter>
bool check_equality(Iter first1, Iter last1, Iter first2) {
  for ( ; first1 != last1; ++first1, ++first2) {
    if (*first1 != *first2) {
      nout << "Difference detected: ";
      nout << "****** ph1 ********* ";
      first1->ascii_dump(nout);
      nout << "****** ph2 ********* ";
      first2->ascii_dump(nout);
      nout << "\n";
      return false;
    }
  }
  return true;
}

bool
computation(dim_type dim, Poly& ph,
            const Cons& hyper_cs, Cons& cons) {
  std::vector<Poly> ph_split;
  std::vector<Poly> knres;
  ph_split.resize(10 + 1, ph);
  knres.resize(10 + 1, ph);

  ph.add_cons(hyper_cs);
  ph.minimize();
  knres[0] = ph;
  ph_split[0] = ph;

  nout << "\nSplitting hypercube of dim " << dim << " | ";
  {
    Clock clock;
    for (dim_type i = 0; i < 10; ++i)
      ph_split[i + 1] = ph_split[i].integral_split(cons[i]);
    clock.print_elapsed(nout);
  }

  nout << "\nUsual eager mode \t     | ";
  {
    Clock clock;
    for (dim_type i = 0; i < 10; ++i) {
      auto [c_then, c_else] = detail::integral_complement_cons(cons[i]);
      knres[i + 1] = knres[i];
      knres[i].add_con(c_then);
      knres[i].minimize();
      knres[i + 1].add_con(c_else);
      knres[i + 1].minimize();
    }
    clock.print_elapsed(nout);
  }

  // Check
  bool ok;
  for (dim_type i = 0; i <= 10; ++i) {
    ok = (ph_split[i] == knres[i]);
    if (!ok) {
      nout << " With constraint : ";
      cons[i].ascii_dump(nout);
      nout << "****** ph_then ********* ";
      ph_split[i].ascii_dump(nout);
      nout << "****** knres_then ********* ";
      knres[i].ascii_dump(nout);
      return ok;
    }
  }
  nout << "\n";
  return ok;
}



bool
test01() {
  dim_type min_dim = 5;
  dim_type max_dim = 8;
  if (check_exp_eval()) {
    min_dim = 10;
    max_dim = 12;
  }

  for (dim_type dim = min_dim; dim <= max_dim; ++dim) {
    Cons hyper_cs;
    hyper_cs.reserve(2 * dim);
    for (dim_type i = 0; i < dim; ++i) {
      hyper_cs.push_back(Var(i) >= 0);
      hyper_cs.push_back(Var(i) <= 10);
    }

    Poly ph(dim, Topol::CLOSED);
    Cons cons;
    cons.reserve(10);
    for (dim_type j = 1; j <= 10; ++j)
      cons.push_back(Var(0) < j);

    if (not computation(dim, ph, hyper_cs, cons))
      return false;
  }
  return true;
}

bool
test02() {
  dim_type min_dim = 5;
  dim_type max_dim = 8;
  if (check_exp_eval()) {
    min_dim = 10;
    max_dim = 12;
  }

  for (dim_type dim = min_dim; dim <= max_dim; ++dim) {
    Cons hyper_cs;
    hyper_cs.reserve(2 * dim);
    for (dim_type i = 0; i < dim; ++i) {
      hyper_cs.push_back(Var(i) >= 0);
      hyper_cs.push_back(Var(i) <= 10);
    }

    Poly ph(dim, Topol::CLOSED);
    Cons cons;
    cons.reserve(10);
    for (dim_type j = 1; j <= 10; ++j)
      cons.push_back(Var(0) <= j);

    if (not computation(dim, ph, hyper_cs, cons))
      return false;
  }

  return true;
}

bool
test03() {
  dim_type min_dim = 5;
  dim_type max_dim = 8;
  if (check_exp_eval()) {
    min_dim = 10;
    max_dim = 12;
  }

  for (dim_type dim = min_dim; dim <= max_dim; ++dim) {
    Cons hyper_cs;
    hyper_cs.reserve(2 * dim);
    for (dim_type i = 0; i < dim; ++i) {
      hyper_cs.push_back(Var(i) >= 0);
      hyper_cs.push_back(Var(i) <= 10);
    }

    Poly ph(dim);
    Cons cons;
    cons.reserve(10);
    for (dim_type j = 1; j <= 10; ++j)
      cons.push_back(Var(0) <= j);

    if (not computation(dim, ph, hyper_cs, cons))
      return false;
  }
  return true;
}

bool
test04() {
  dim_type dim = 5;
  dim_type split_vars = 3;
  if (check_exp_eval()) {
    dim = 8;
    split_vars = 8;
  }

  Cons hyper_cs;
  hyper_cs.reserve(2 * dim);
  for (dim_type i = 0; i < dim; ++i) {
    hyper_cs.push_back(Var(i) >= 0);
    hyper_cs.push_back(Var(i) <= 8);
  }

  Poly ph_init(dim, Topol::CLOSED);
  ph_init.add_cons(hyper_cs);
  ph_init.minimize();

  const dim_type num_poly = 1 << split_vars;
  auto phs_old = std::vector<Poly>(num_poly, Poly(0, Topol::CLOSED));
  auto phs_new = phs_old;

  nout << "New splitting of hypercube | ";
  {
    auto& phs = phs_new;
    Clock clock;
    phs[0] = ph_init;
    for (dim_type sv = 0; sv < split_vars; ++sv) {
      Con sc { Var(sv) <= 4 };
      const dim_type delta = 1 << sv;
      for (dim_type i = 0; i < delta; ++i)
        phs[i + delta] = phs[i].integral_split(sc);
    }
    for (auto& ph : phs)
      ph.minimize();
    clock.print_elapsed(nout);
    nout << std::endl;
  }

  nout << "Old splitting of hypercube | ";
  {
    auto& phs = phs_old;
    Clock clock;
    phs[0] = ph_init;
    for (dim_type sv = 0; sv < split_vars; ++sv) {
      Con sc { Var(sv) <= 4 };
      auto [c_then, c_else] = detail::integral_complement_cons(sc);
      const dim_type delta = 1 << sv;
      for (dim_type i = 0; i < delta; ++i) {
        phs[i + delta] = phs[i];
        phs[i].add_con(c_then);
        phs[i + delta].add_con(c_else);
        // Being eager is actually better.
        phs[i].minimize();
        phs[i + delta].minimize();
      }
    }
    for (auto& ph : phs)
      ph.minimize();
    clock.print_elapsed(nout);
    nout << std::endl;
  }

  return check_equality(phs_old.begin(), phs_old.end(), phs_new.begin());
}

bool
test05() {
  dim_type min_dim = 5;
  dim_type max_dim = 8;
  if (check_exp_eval()) {
    min_dim = 8;
    max_dim = 12;
  }

  for (dim_type dim = min_dim; dim <= max_dim; ++dim) {
    Cons hyper_cs;
    hyper_cs.reserve(2 * dim);
    for (dim_type i = 0; i < dim; ++i) {
      hyper_cs.push_back(Var(i) >= 0);
      hyper_cs.push_back(Var(i) <= 10);
    }

    Poly ph(dim, Topol::CLOSED);
    Cons cons;
    cons.reserve(10);
    for (dim_type j = 1; j <= 10; ++j)
      cons.push_back(Var(0) <= j);

    if (not computation(dim, ph, hyper_cs, cons))
      return false;
  }
  return true;
}

} // namespace

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
  DO_TEST(test03);
  DO_TEST(test04);
  DO_TEST(test05);
END_MAIN
