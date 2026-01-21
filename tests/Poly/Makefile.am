# PPLite: a lightweight library for convex polyhedra derived from PPL.
# Copyright (C) 2001-2010 Roberto Bagnara <bagnara@cs.unipr.it>
# Copyright (C) 2010-2018 BUGSENG srl (http://bugseng.com)
# Copyright (C) 2018-2026 Enea Zaffanella <enea.zaffanella@unipr.it>
#
# This file is part of PPLite.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

include $(top_srcdir)/tests/common_preamble.am

XFAIL_TESTS =

# List tests here.

if THREAD_SAFE
THREAD_SAFETY_TESTS = threadsafe
else !THREAD_SAFE
THREAD_SAFETY_TESTS =
endif !THREAD_SAFE

TESTS = \
$(THREAD_SAFETY_TESTS) \
addcons1 \
addcons2 \
addcons3 \
addgens1 \
addgens2 \
addspacedims \
assign_ops \
affinedim1 \
affineimage \
affinepreimage \
bbox1 \
bbox2 \
boxed_inters \
bounded \
concatenate \
conhull \
constrains \
disjoint \
dualhyper \
efc \
empty1 \
equals \
expandspacedim \
foldspacedims \
gentocon \
hash \
hyper \
hyperphull \
inters \
mapspacedims \
mater \
minmax \
parallelimage \
polydiff \
polyhull \
randphull \
relwith1 \
relwith2 \
removehigherdims \
removespacedims \
split \
split2 \
split_dualhyper \
split_hyper \
split_hyper2 \
split_int \
split_int_dualhyper \
split_int_hyper \
split_int_rand \
split_rand \
timeelapse \
topclosure \
univ1 \
widen_bhrz03_1 \
widen_bhrz03_2 \
widen_bhrz03_3 \
widen_h79 \
widen_sas18

# List tests' sources here.

addcons1_SOURCES = addcons1.cc
addcons2_SOURCES = addcons2.cc
addcons3_SOURCES = addcons3.cc
addgens1_SOURCES = addgens1.cc
addgens2_SOURCES = addgens2.cc
addspacedims_SOURCES = addspacedims.cc
assign_ops_SOURCES = assign_ops.cc
affinedim1_SOURCES = affinedim1.cc
affineimage_SOURCES = affineimage.cc
affinepreimage_SOURCES = affinepreimage.cc
bbox1_SOURCES = bbox1.cc
bbox2_SOURCES = bbox2.cc
boxed_inters_SOURCES = boxed_inters.cc
bounded_SOURCES = bounded.cc
concatenate_SOURCES = concatenate.cc
conhull_SOURCES = conhull.cc
constrains_SOURCES = constrains.cc
disjoint_SOURCES = disjoint.cc
dualhyper_SOURCES = dualhyper.cc
efc_SOURCES = efc.cc
empty1_SOURCES = empty1.cc
equals_SOURCES = equals.cc
expandspacedim_SOURCES = expandspacedim.cc
foldspacedims_SOURCES = foldspacedims.cc
gentocon_SOURCES = gentocon.cc
hash_SOURCES = hash.cc
hyper_SOURCES = hyper.cc
hyperphull_SOURCES = hyperphull.cc
inters_SOURCES = inters.cc
mapspacedims_SOURCES = mapspacedims.cc
mater_SOURCES = mater.cc
minmax_SOURCES = minmax.cc
parallelimage_SOURCES = parallelimage.cc
polydiff_SOURCES = polydiff.cc
polyhull_SOURCES = polyhull.cc
randphull_SOURCES = randphull.cc
relwith1_SOURCES = relwith1.cc
relwith2_SOURCES = relwith2.cc
removehigherdims_SOURCES = removehigherdims.cc
removespacedims_SOURCES = removespacedims.cc
split_SOURCES = split.cc
split2_SOURCES = split2.cc
split_dualhyper_SOURCES = split_dualhyper.cc
split_hyper_SOURCES = split_hyper.cc
split_hyper2_SOURCES = split_hyper2.cc
split_int_SOURCES = split_int.cc
split_int_dualhyper_SOURCES = split_int_dualhyper.cc
split_int_hyper_SOURCES = split_int_hyper.cc
split_int_rand_SOURCES = split_int_rand.cc
split_rand_SOURCES = split_rand.cc
threadsafe_SOURCES = threadsafe.cc
timeelapse_SOURCES = timeelapse.cc
topclosure_SOURCES = topclosure.cc
univ1_SOURCES = univ1.cc
widen_bhrz03_1_SOURCES = widen_bhrz03_1.cc
widen_bhrz03_2_SOURCES = widen_bhrz03_2.cc
widen_bhrz03_3_SOURCES = widen_bhrz03_3.cc
widen_h79_SOURCES = widen_h79.cc
widen_sas18_SOURCES = widen_sas18.cc


# Checking Dyn_Poly kinds (using the same tests)

DYN_POLY_KINDS = \
POLY \
B_POLY \
U_POLY \
F_POLY \
UF_POLY

# default value
KIND=F_POLY

.PHONY: check_all_dyn_poly

check_all_dyn_poly:
	+@for kind in ${DYN_POLY_KINDS}; do \
	  make -s clean; \
	  echo "====="; \
	  echo "===== Running Dyn_Poly tests for kind = $$kind"; \
	  echo "====="; \
	  make V=0 check TEST_CPPFLAGS=-DTEST_ABS_POLY_KIND=$$kind || break; \
	  echo "====="; \
	  echo "===== Completed Dyn_Poly tests for kind = $$kind" ; \
	  echo "====="; \
	done

check_dyn_poly:
	+@echo -e "\n===== Admissible values for KIND variable:"; \
	echo -e "===== ${DYN_POLY_KINDS}\n"; \
	make -s clean; \
	echo "====="; \
	echo "===== Running Dyn_Poly tests for kind = ${KIND}"; \
	echo "====="; \
	make V=1 check TEST_CPPFLAGS=-DTEST_ABS_POLY_KIND=${KIND}; \
	echo "====="; \
	echo "===== Completed Dyn_Poly tests for kind = ${KIND}" ; \
	echo "=====";
