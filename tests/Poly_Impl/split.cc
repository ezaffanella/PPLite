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

bool
test01() {
  Var A(0);
  Var B(1);
  Var C(2);
  Var D(3);

  Poly ph;
  std::string input = R"###(topol NNC
dim 4
status MINIMIZED
=> cs sys
sing_rows 0
sk_rows 16
> : dim 4 : -1 -1 1 1  : 4
> : dim 4 : -1 1 -1 1  : 4
> : dim 4 : -1 1 1 1  : 0
> : dim 4 : -1 1 1 -1  : 4
> : dim 4 : 1 -1 -1 1  : 4
> : dim 4 : 1 -1 1 1  : 0
> : dim 4 : 1 -1 1 -1  : 4
> : dim 4 : 1 1 -1 1  : 0
>= : dim 4 : 1 1 -1 -1  : 4
> : dim 4 : 1 1 1 1  : -4
> : dim 4 : 1 1 1 -1  : 0
>= : dim 4 : -1 -1 -1 -1  : 9
> : dim 4 : -1 -1 -1 1  : 7
> : dim 4 : -1 -1 1 -1  : 7
>= : dim 4 : -1 1 -1 -1  : 7
>= : dim 4 : 1 -1 -1 -1  : 7
ns_rows 0
=> gs sys
sing_rows 0
sk_rows 36
C : dim 4 : 4 4 3 -3  : 2
C : dim 4 : 4 3 4 -3  : 2
C : dim 4 : 3 4 4 -3  : 2
C : dim 4 : 4 3 9 2  : 2
C : dim 4 : 3 4 9 2  : 2
C : dim 4 : 4 9 3 2  : 2
C : dim 4 : 3 9 4 2  : 2
C : dim 4 : 9 4 3 2  : 2
C : dim 4 : 9 3 4 2  : 2
C : dim 4 : 4 3 -3 4  : 2
C : dim 4 : 3 4 -3 4  : 2
C : dim 4 : 4 3 2 9  : 2
C : dim 4 : 3 4 2 9  : 2
C : dim 4 : 3 9 2 4  : 2
C : dim 4 : 9 3 2 4  : 2
C : dim 4 : 4 4 -3 3  : 2
C : dim 4 : 4 9 2 3  : 2
C : dim 4 : 9 4 2 3  : 2
C : dim 4 : 3 -3 4 4  : 2
C : dim 4 : 4 -3 4 3  : 2
C : dim 4 : 4 -3 3 4  : 2
C : dim 4 : 3 2 4 9  : 2
C : dim 4 : 4 2 3 9  : 2
C : dim 4 : 3 2 9 4  : 2
C : dim 4 : 4 2 9 3  : 2
C : dim 4 : 9 2 4 3  : 2
C : dim 4 : 9 2 3 4  : 2
C : dim 4 : -3 4 4 3  : 2
C : dim 4 : -3 4 3 4  : 2
C : dim 4 : -3 3 4 4  : 2
C : dim 4 : 2 4 3 9  : 2
C : dim 4 : 2 3 4 9  : 2
C : dim 4 : 2 4 9 3  : 2
C : dim 4 : 2 3 9 4  : 2
C : dim 4 : 2 9 4 3  : 2
C : dim 4 : 2 9 3 4  : 2
ns_rows 2
2 : { 21, 23 }
2 : { 31, 33 }
sat_c
36 x 16
0 1 0 1 1 0 1 1 1 0 1 1 0 1 1 1 
1 0 0 1 1 1 1 0 1 0 1 1 0 1 1 1 
1 1 1 1 0 0 1 0 1 0 1 1 0 1 1 1 
1 0 1 1 1 1 1 0 1 1 1 0 0 1 1 1 
1 1 1 1 0 1 1 0 1 1 1 0 0 1 1 1 
0 1 1 1 1 0 1 1 1 1 1 0 0 1 1 1 
1 1 1 1 0 0 1 1 1 1 1 0 0 1 1 1 
0 1 0 1 1 1 1 1 1 1 1 0 0 1 1 1 
1 0 0 1 1 1 1 1 1 1 1 0 0 1 1 1 
1 1 0 0 1 1 1 1 1 0 0 1 1 0 1 1 
1 1 1 1 1 0 0 1 1 0 0 1 1 0 1 1 
1 1 1 0 1 1 1 1 1 1 0 0 1 0 1 1 
1 1 1 1 1 1 0 1 1 1 0 0 1 0 1 1 
1 1 1 1 1 0 0 1 1 1 1 0 1 0 1 1 
1 1 0 0 1 1 1 1 1 1 1 0 1 0 1 1 
0 1 0 1 1 0 1 1 1 0 1 1 1 0 1 1 
0 1 1 1 1 0 1 1 1 1 1 0 1 0 1 1 
0 1 0 1 1 1 1 1 1 1 1 0 1 0 1 1 
1 1 1 1 1 1 1 0 0 0 0 1 1 1 0 1 
1 0 0 1 1 1 1 0 1 0 1 1 1 1 0 1 
1 1 0 0 1 1 1 1 1 0 0 1 1 1 0 1 
1 1 1 1 1 1 1 1 0 1 0 0 1 1 0 1 
1 1 1 0 1 1 1 1 1 1 0 0 1 1 0 1 
1 1 1 1 1 1 1 0 0 1 1 0 1 1 0 1 
1 0 1 1 1 1 1 0 1 1 1 0 1 1 0 1 
1 0 0 1 1 1 1 1 1 1 1 0 1 1 0 1 
1 1 0 0 1 1 1 1 1 1 1 0 1 1 0 1 
1 1 1 1 0 0 1 0 1 0 1 1 1 1 1 0 
1 1 1 1 1 0 0 1 1 0 0 1 1 1 1 0 
1 1 1 1 1 1 1 0 0 0 0 1 1 1 1 0 
1 1 1 1 1 1 0 1 1 1 0 0 1 1 1 0 
1 1 1 1 1 1 1 1 0 1 0 0 1 1 1 0 
1 1 1 1 0 1 1 0 1 1 1 0 1 1 1 0 
1 1 1 1 1 1 1 0 0 1 1 0 1 1 1 0 
1 1 1 1 0 0 1 1 1 1 1 0 1 1 1 0 
1 1 1 1 1 0 0 1 1 1 1 0 1 1 1 0 
sat_g
16 x 36
0 1 1 1 1 0 1 0 1 1 1 1 1 1 1 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
1 0 1 0 1 1 1 1 0 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 0 0 1 1 1 1 1 1 1 1 1 1 
0 0 1 1 1 1 1 0 0 0 1 1 1 1 0 0 1 0 1 0 0 1 1 1 1 0 0 1 1 1 1 1 1 1 1 1 
1 1 1 1 1 1 1 1 1 0 1 0 1 1 0 1 1 1 1 1 0 1 0 1 1 1 0 1 1 1 1 1 1 1 1 1 
1 1 0 1 0 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 0 1 0 1 
0 1 0 1 1 0 0 1 1 1 0 1 1 0 1 0 0 1 1 1 1 1 1 1 1 1 1 0 0 1 1 1 1 1 0 0 
1 1 1 1 1 1 1 1 1 1 0 1 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 0 1 1 1 1 0 
1 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 1 1 1 0 0 1 1 0 1 0 1 1 0 0 1 1 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 0 1 0 1 1 1 1 1 0 1 0 1 0 1 1 
0 0 0 1 1 1 1 1 1 0 0 1 1 1 1 0 1 1 0 0 0 1 1 1 1 1 1 0 0 0 1 1 1 1 1 1 
1 1 1 1 1 1 1 1 1 0 0 0 0 1 1 1 1 1 0 1 0 0 0 1 1 1 1 1 0 0 0 0 1 1 1 1 
1 1 1 0 0 0 0 0 0 1 1 0 0 0 0 1 0 0 1 1 1 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 
=> cs_pending
sing_rows 0
sk_rows 0
ns_rows 0
=> gs_pending
sing_rows 0
sk_rows 0
ns_rows 0
)###";
  std::stringstream ss_in(input);
  ph.ascii_load(ss_in);

  Con c(B + C + D < 8);

  Poly knres_then(ph);
  knres_then.add_con(c);
  Poly ph_else = ph.split(c);

  bool ok = (knres_then == ph);
  return ok;
}

BEGIN_MAIN
  DO_TEST(test01);
END_MAIN

