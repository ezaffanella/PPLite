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

#include "pplite_test.hh"

bool
test01() {
  Poly ph;
  std::stringstream ss;
  ph.ascii_dump(ss);
  auto expected = R"###(topol C
dim 0
status MINIMIZED
=> cs sys
sing_rows 0
sk_rows 1
> : dim 0 :  : 1
ns_rows 0
=> gs sys
sing_rows 0
sk_rows 1
P : dim 0 :  : 1
ns_rows 0
sat_c
1 x 1
1 
sat_g
1 x 1
1 
=> cs_pending
sing_rows 0
sk_rows 0
ns_rows 0
=> gs_pending
sing_rows 0
sk_rows 0
ns_rows 0
)###";
  return check_print(ss.str(), expected);
}

bool
test02() {
  Poly ph(3);
  std::stringstream ss;
  ph.ascii_dump(ss);
  auto expected = R"###(topol C
dim 3
status MINIMIZED
=> cs sys
sing_rows 0
sk_rows 1
> : dim 0 :  : 1
ns_rows 0
=> gs sys
sing_rows 3
L : dim 1 : 1  : 0
L : dim 2 : 0 1  : 0
L : dim 3 : 0 0 1  : 0
sk_rows 1
P : dim 0 :  : 1
ns_rows 0
sat_c
1 x 1
1 
sat_g
1 x 1
1 
=> cs_pending
sing_rows 0
sk_rows 0
ns_rows 0
=> gs_pending
sing_rows 0
sk_rows 0
ns_rows 0
)###";
  return check_print(ss.str(), expected);
}

bool
test03() {
  Poly ph;
  std::string input = R"###(topol NNC
dim 3
status PENDING
=> cs sys
sing_rows 1
= : dim 1 : 1  : 0
sk_rows 2
> : dim 2 : 0 -562949953421312  : -3537115888337719
> : dim 2 : 0 0  : 1
ns_rows 0
=> gs sys
sing_rows 1
L : dim 3 : 0 0 1  : 0
sk_rows 2
C : dim 2 : 0 -3537115888337719  : 562949953421312
R : dim 2 : 0 -1  : 0
ns_rows 1
2 : { 0, 1 }
sat_c
2 x 2
0 1 
1 0 
sat_g
2 x 2
0 1 
1 0 
=> cs_pending
sing_rows 0
sk_rows 1
> : dim 2 : 0 -281474976710656  : -3537115888337719
ns_rows 0
=> gs_pending
sing_rows 0
sk_rows 0
ns_rows 0
)###";

  std::stringstream ss_in(input);
  ph.ascii_load(ss_in);
  std::stringstream ss_out;
  ph.ascii_dump(ss_out);
  if (input != ss_out.str())
    return false;
  if (!ph.check_inv())
    return false;
  ph.minimize();
  return ph.check_inv();
}

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
  DO_TEST(test03);
END_MAIN

