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

bool test01() {
  Var A(0), B(1), C(2), D(3), E(4), F(5), G(6), H(7), I(8), J(9),
    K(10), L(11), M(12), N(13), O(14), P(15), Q(16);

  Cons cs1 = {
    Q == 0, P == 0, O == 0, N == 0, M == 0, H == 0, G == 0, F == 0, E == 0,
    C - D + K == 0, B - C + J == 0, -A >= -1, -8*B + 8*C >= 5,
    -B >= -4, -16*C + 16*D >= 15, -C >= -8, -D >= -12,
    -L >= -16, 4*L >= 5, I >= 0, B - I >= 0, A >= 0
  };

  Cons cs2 = {
    Q == 0, P == 0, O == 0, N == 0, M == 0, H == 0, G == 0, F == 0, E == 0,
    B == 0, A == 1, -C >= -8, -C + D >= 0, -D >= -12, -I >= -4, -J >= -8,
    -K >= -12, -L >= -16, 4*L >= 5, 16*K >= 15, 8*J >= 5, I >= 0, C >= 0
  };

  const auto sd = 17;
  Poly ph1(sd, Topol::NNC);
  ph1.add_cons(cs1);
  ph1.minimize();

  Poly ph2(sd, Topol::NNC);
  ph2.add_cons(cs2);
  ph2.minimize();

  std::string input1 = R"###(
dim 17
empty 0
topol NNC
is_normalized 0
===itvs-component===
0 : block-dim
1 : block-dim
2 : block-dim
3 : block-dim
4 : [ lb 0 , ub 0 ]
5 : [ lb 0 , ub 0 ]
6 : [ lb 0 , ub 0 ]
7 : [ lb 0 , ub 0 ]
8 : block-dim
9 : block-dim
10 : block-dim
11 : [ lb 5/4 , ub 16 ]
12 : [ lb 0 , ub 0 ]
13 : [ lb 0 , ub 0 ]
14 : [ lb 0 , ub 0 ]
15 : [ lb 0 , ub 0 ]
16 : [ lb 0 , ub 0 ]
===blocks-component===
blocks 1
7 : { 8 9 2 3 10 1 0 }
factors 1
===start-of-factor 0
( block = 7 : { 8 9 2 3 10 1 0 } )
topol NNC
dim 7
status MINIMIZED
=> cs sys
sing_rows 2
= : dim 7 : 0 1 -1 0 0 1 0  : 0
= : dim 7 : 0 0 1 -1 1 0 0  : 0
sk_rows 9
>= : dim 7 : -1 -1 1 0 0 0 0  : 0
>= : dim 7 : 0 0 0 0 0 0 1  : 0
>= : dim 7 : 0 0 0 0 0 0 -1  : 1
>= : dim 5 : 0 8 0 0 0  : -5
>= : dim 1 : 1  : 0
>= : dim 6 : 0 0 0 -1 0 0  : 12
>= : dim 7 : 0 0 -16 16 0 0 0  : -15
>= : dim 7 : 0 1 -1 0 0 0 0  : 4
>= : dim 7 : 0 0 -1 0 0 0 0  : 8
ns_rows 1
9 : { 0, 1, 2, 3, 4, 5, 6, 7, 8 }
=> gs sys
sing_rows 0
sk_rows 24
P : dim 5 : 0 10 10 25 15  : 16
P : dim 7 : 0 10 10 25 15 0 16  : 16
P : dim 6 : 64 10 74 89 15 64  : 16
P : dim 7 : 64 10 74 89 15 64 16  : 16
P : dim 5 : 0 128 128 143 15  : 16
P : dim 7 : 0 128 128 143 15 0 16  : 16
P : dim 6 : 64 64 128 143 15 64  : 16
P : dim 7 : 64 64 128 143 15 64 16  : 16
P : dim 6 : 0 5 37 96 59 32  : 8
P : dim 5 : 0 5 5 96 91  : 8
P : dim 7 : 0 5 37 96 59 32 8  : 8
P : dim 7 : 0 5 5 96 91 0 8  : 8
P : dim 6 : 32 5 37 96 59 32  : 8
P : dim 7 : 32 5 37 96 59 32 8  : 8
P : dim 5 : 0 8 8 12 4  : 1
P : dim 7 : 0 8 8 12 4 0 1  : 1
P : dim 6 : 4 4 8 12 4 4  : 1
P : dim 7 : 4 4 8 12 4 4 1  : 1
P : dim 6 : 0 10 74 89 15 64  : 16
P : dim 7 : 0 10 74 89 15 64 16  : 16
P : dim 6 : 0 4 8 12 4 4  : 1
P : dim 7 : 0 4 8 12 4 4 1  : 1
P : dim 6 : 0 64 128 143 15 64  : 16
P : dim 7 : 0 64 128 143 15 64 16  : 16
ns_rows 0
sat_c
24 x 9
0 0 1 0 0 1 0 1 1 
0 1 0 0 0 1 0 1 1 
0 0 1 0 1 1 0 0 1 
0 1 0 0 1 1 0 0 1 
0 0 1 1 0 1 0 1 0 
0 1 0 1 0 1 0 1 0 
0 0 1 1 1 1 0 0 0 
0 1 0 1 1 1 0 0 0 
1 0 1 0 0 0 1 0 1 
0 0 1 0 0 0 1 1 1 
1 1 0 0 0 0 1 0 1 
0 1 0 0 0 0 1 1 1 
0 0 1 0 1 0 1 0 1 
0 1 0 0 1 0 1 0 1 
0 0 1 1 0 0 1 1 0 
0 1 0 1 0 0 1 1 0 
0 0 1 1 1 0 1 0 0 
0 1 0 1 1 0 1 0 0 
1 0 1 0 0 1 0 0 1 
1 1 0 0 0 1 0 0 1 
1 0 1 1 0 0 1 0 0 
1 1 0 1 0 0 1 0 0 
1 0 1 1 0 1 0 0 0 
1 1 0 1 0 1 0 0 0 
sat_g
9 x 24
0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 1 1 1 1 1 1 
0 1 0 1 0 1 0 1 0 0 1 1 0 1 0 1 0 1 0 1 0 1 0 1 
1 0 1 0 1 0 1 0 1 1 0 0 1 0 1 0 1 0 1 0 1 0 1 0 
0 0 0 0 1 1 1 1 0 0 0 0 0 0 1 1 1 1 0 0 1 1 1 1 
0 0 1 1 0 0 1 1 0 0 0 0 1 1 0 0 1 1 0 0 0 0 0 0 
1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 1 1 0 0 1 1 
0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 0 0 1 1 0 0 
1 1 0 0 1 1 0 0 0 1 0 1 0 0 1 1 0 0 0 0 0 0 0 0 
1 1 1 1 0 0 0 0 1 1 1 1 1 1 0 0 0 0 1 1 0 0 0 0 
=> cs_pending
sing_rows 0
sk_rows 0
ns_rows 0
=> gs_pending
sing_rows 0
sk_rows 0
ns_rows 0
===end-of-factor 0
)###";

  std::string input2 = R"###(
dim 17
empty 0
topol NNC
is_normalized 0
===itvs-component===
0 : [ lb 1 , ub 1 ]
1 : [ lb 0 , ub 0 ]
2 : block-dim
3 : block-dim
4 : [ lb 0 , ub 0 ]
5 : [ lb 0 , ub 0 ]
6 : [ lb 0 , ub 0 ]
7 : [ lb 0 , ub 0 ]
8 : [ lb 0 , ub 4 ]
9 : block-dim
10 : [ lb 15/16 , ub 12 ]
11 : [ lb 5/4 , ub 16 ]
12 : [ lb 0 , ub 0 ]
13 : [ lb 0 , ub 0 ]
14 : [ lb 0 , ub 0 ]
15 : [ lb 0 , ub 0 ]
16 : [ lb 0 , ub 0 ]
===blocks-component===
blocks 1
3 : { 2 3 9 }
factors 1
===start-of-factor 0
( block = 3 : { 2 3 9 } )
topol NNC
dim 3
status MINIMIZED
=> cs sys
sing_rows 0
sk_rows 6
>= : dim 3 : -1 1 0  : 0
>= : dim 3 : 1 0 0  : 0
>= : dim 3 : -1 0 0  : 8
>= : dim 3 : 0 -1 0  : 12
>= : dim 3 : 0 0 -1  : 8
>= : dim 3 : 0 0 8  : -5
ns_rows 1
6 : { 0, 1, 2, 3, 4, 5 }
=> gs sys
sing_rows 0
sk_rows 8
P : dim 3 : 0 0 8  : 1
P : dim 3 : 0 12 8  : 1
P : dim 3 : 8 8 8  : 1
P : dim 3 : 8 12 8  : 1
P : dim 3 : 0 0 5  : 8
P : dim 3 : 0 96 5  : 8
P : dim 3 : 64 64 5  : 8
P : dim 3 : 64 96 5  : 8
ns_rows 0
sat_c
8 x 6
0 0 1 1 0 1 
1 0 1 0 0 1 
0 1 0 1 0 1 
1 1 0 0 0 1 
0 0 1 1 1 0 
1 0 1 0 1 0 
0 1 0 1 1 0 
1 1 0 0 1 0 
sat_g
6 x 8
0 1 0 1 0 1 0 1 
0 0 1 1 0 0 1 1 
1 1 0 0 1 1 0 0 
1 0 1 0 1 0 1 0 
0 0 0 0 1 1 1 1 
1 1 1 1 0 0 0 0 
=> cs_pending
sing_rows 0
sk_rows 0
ns_rows 0
=> gs_pending
sing_rows 0
sk_rows 0
ns_rows 0
===end-of-factor 0
)###";

  std::stringstream ss1(input1);
  F_Poly fph1;
  if (not fph1.ascii_load(ss1))
    return false;

  std::stringstream ss2(input2);
  F_Poly fph2;
  if (not fph2.ascii_load(ss2))
    return false;

  Vars vars = { A, E, F, G, H, I , K, L, M, N, O, P, Q };
  Index_Set iset;
  for (auto v : vars)
    iset.set(v.id());

  ph1.remove_space_dims(iset);
  ph2.remove_space_dims(iset);
  fph1.remove_space_dims(iset);
  fph2.remove_space_dims(iset);

  ph1.minimize();
  ph2.minimize();
  fph1.minimize();
  fph2.minimize();

#if 1
  // This permutation makes things work.
  Dims perm = { 3, 1, 2, 0 };
  ph1.map_space_dims(perm);
  ph2.map_space_dims(perm);
  fph1.map_space_dims(perm);
  fph2.map_space_dims(perm);
#endif

  nout << "\n===== ph1 =====\n";
  ph1.print(nout);
  nout << "\n===== fph1 =====\n";
  fph1.print(nout);
  nout << "\n=====\n";

  nout << "\n===== ph2 =====\n";
  ph2.print(nout);
  nout << "\n===== fph2 =====\n";
  fph2.print(nout);
  nout << "\n=====\n";

  Poly ph3 = ph1;
  ph3.con_hull_assign(ph2);
  ph3.minimize();

  F_Poly fph3 = fph1;
  fph3.con_hull_assign(fph2);
  fph3.minimize();

  nout << "\n===== ph3 =====\n";
  ph3.print(nout);
  nout << "\n===== fph3 =====\n";
  fph3.print(nout);

  if (ph1 != fph1.to_poly()) {
    nout << "\nph1 and fph1 differ\n";
    return false;
  }
  if (ph2 != fph2.to_poly()) {
    nout << "\nph2 and fph2 differ\n";
    return false;
  }
  if (ph3 != fph3.to_poly()) {
    nout << "\nout3 and fph3 differ\n";
    return false;
  }
  return true;
}

bool
test02() {
  Var A(0); Var B(1); Var C(2); Var D(3); Var E(4);
  Var F(5); Var G(6); Var H(7); Var I(8); Var J(9);
  Var K(10); Var L(11); Var M(12); Var N(13); Var O(14);
  Var P(15); Var Q(16);

  Cons cs_x = { A == 1, B == 1, C + K - L == 12, J + K == 12, E == 1, F == 2,
                G == 0, H == 0, I == 0, M == 0, N == 0, O == 0, P == 0, Q == 0,
                D == 12, -K >= -12, C + K >= 12, -C >= -8 };
  F_Poly x(17);
  x.add_cons(std::move(cs_x));

  Cons cs_y = { A == 1, D == 12, L == 0, E == 1, F == 2, G == 0, H == 0, M == 0,
                N == 0, O == 0, P == 0, Q == 0, B >= 1, -B >= -4, C >= 0, J >= 0,
                -C >= -8, -J >= -8, -K >= -12, 16*K >= 15, I >= 0, -I >= -4 };
  F_Poly y(17);
  y.add_cons(std::move(cs_y));

  auto kr = x;
  kr.con_hull_assign(y);
  kr.minimize();

  con_hull(x, &y, &y + 1);
  x.minimize();

  bool ok = (kr == x);

  return ok;
}

bool
test03() {
  Var A(0); Var B(1); Var C(2); Var D(3);

  Cons cs_x = { A + C - D == 12, B + C == 12,
                -C >= -12, A + C >= 12, -A >= -8 };
  F_Poly x(4);
  x.add_cons(std::move(cs_x));

  Cons cs_y = {  D == 0, A >= 0, B >= 0,
                -A >= -8, -B >= -8, -C >= -12, 16*C >= 15 };
  F_Poly y(4);
  y.add_cons(std::move(cs_y));

  auto kr = x;
  kr.con_hull_assign(y);

  auto nary = x;
  con_hull(nary, &y, &y + 1);

  bool ok = (kr == nary);

  return ok;
}

bool
test04() {
  Var A(0); Var B(1); Var C(2); Var D(3);

  Cons cs_x = { A - C - D == 12, B + C == 12,
                -C >= -12, A >= 0, -A >= -8  };
  F_Poly fx(4);
  fx.add_cons(cs_x);
  fx.minimize();
  Poly x(4);
  x.add_cons(std::move(cs_x));
  x.minimize();

  Cons cs_y = {  D == 0, A >= 0, B >= 0, -A >= -8, -B >= -8,
                 -C >= -12, 16*C >= 15 };
  F_Poly fy(4);
  fy.add_cons(cs_y);
  fy.minimize();
  Poly y(4);
  y.add_cons(std::move(cs_y));
  y.minimize();

  auto fchull = fx;
  con_hull(fchull, &fy, &fy + 1);
  fchull.minimize();
  auto chull = x;
  con_hull(chull, &y, &y + 1);
  chull.minimize();

  bool ok = (chull == fchull.to_poly());

  return ok;
}

#if 0 // TEST for the old specialized con_hull implementation
// Showing that it yields a different incomparable
// result than the same operator on the corresponding Polys,
// due to the missing normalization of the factorization.
bool
test05() {

  Var A(0); Var B(1); Var C(2); Var D(3); Var E(4);
  Var F(5); Var G(6); Var H(7); Var I(8); Var J(9);
  Var K(10); Var L(11); Var M(12); Var N(13); Var O(14);
  Var P(15); Var Q(16);

  Cons cs_x = { E == 0, F == 0, G == 0, H == 0, M == 0, N == 0, O == 0, P == 0, Q == 0,
                A >= 0, -A >= -1, B >= 0, -B >= -4, -C + D >= 0, -C >= -8,
                -D >= -12, -K >= -12, 16*K >= 15, -J >= -8, 8*J >= 5,
                -B + C >= 0, I >= 0, -I >= -4, -L >= -16, 4*L >= 5 };

  F_Poly fx(17);
  fx.add_cons(cs_x);
  fx.minimize();

  F_Poly fx1 = fx;

  Poly x(17);
  x.add_cons(std::move(cs_x));
  x.minimize();

  Cons cs_y = { A == 1, B + J - K == 8, I + J == 8, C == 0, E == 0, F == 0, G == 0,
                H == 0, M == 0, N == 0, O == 0, P == 0, Q == 0, -16*B + 16*D - 16*J >= 15,
                -J >= -8, -B >= -4, -D >= -12, 16*B + 16*J >= 143, -L >= -16, 4*L >= 5 };
  F_Poly fy(17);
  fy.add_cons(cs_y);
  fy.minimize();
  F_Poly fy1 = fy;

  Poly y(17);
  y.add_cons(std::move(cs_y));
  y.minimize();

  // lub-block-wise normalization
  auto fchull = fx;
  fchull.con_hull_assign(fy);
  fchull.minimize();

  // Global normalization
  auto fchull1 = fx1;
  con_hull(fchull1, &fy1, &fy1 + 1);
  fchull1.minimize();

  // Poly result
  auto chull = x;
  chull.con_hull_assign(y);
  chull.minimize();

  bool ok = (chull == fchull1.to_poly()
             // local normalization obtains a different result.
             && chull != fchull.to_poly()
             // results are incomparable.
             && !fchull.contains(fchull1)
             && !fchull1.contains(fchull));

  return ok;
}
#endif

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
  DO_TEST(test03);
  DO_TEST(test04);
END_MAIN
