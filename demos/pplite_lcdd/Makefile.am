# PPLite: a lightweight library for convex polyhedra derived from PPL.
# Copyright (C) 2001-2010 Roberto Bagnara <bagnara@cs.unipr.it>
# Copyright (C) 2010-2018 BUGSENG srl (http://bugseng.com)
# Copyright (C) 2018-2023 Enea Zaffanella <enea.zaffanella@unipr.it>
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

SUBDIRS = examples

AM_CPPFLAGS = \
-I$(top_builddir)/src \
-I$(top_srcdir)/src \
@extra_includes@

bin_PROGRAMS = pplite_lcdd

pplite_lcdd_SOURCES = \
pplite_lcdd.cc

pplite_lcdd_LDADD = \
$(top_builddir)/src/libpplite.la \
@extra_libraries@

EXTRA_DIST = \
expected_mpz \
expected_mpz_a

if HAVE_MD5SUM

if ASSERTIONS_ENABLED

TOUGH_EXAMPLES_IF_ASSERTIONS_ENABLED = \
ccc6.ext \
ccp6.ext \
cut32_16.ext \
cyclic16-10.ext \
reg600-5_m.ext \
cyclic17_8.ine \
in6.ine \
in7.ine \
kq20_11_m.ine \
mit31-20.ine \
sampleh8.ine \
trunc10.ine

EXPECTED = expected_mpz_a

else !ASSERTIONS_ENABLED

EXPECTED = expected_mpz

endif !ASSERTIONS_ENABLED

TOUGH_EXAMPLES = \
cp7.ext \
cyclic25_13.ext \
prodst62.ext \
cp6.ine \
mit.ine \
mit288-281.ine \
mit41-16.ine \
mit708-9.ine \
mit71-61.ine \
mit90-86.ine \
$(TOUGH_EXAMPLES_IF_ASSERTIONS_ENABLED)

LCDD = pplite_lcdd$(EXEEXT)

check-local: $(LCDD)
	echo $(TOUGH_EXAMPLES) >tough_examples
	echo $(ECHO_N) "$(ECHO_C)" >obtained
	export LC_ALL=C; \
	for path in $(srcdir)/examples/*.ext $(srcdir)/examples/*.ine; \
	do \
	  file=`basename $$path`; \
	  $(FGREP) $$file tough_examples >/dev/null 2>&1 && continue; \
	  echo $(ECHO_N) $$file " $(ECHO_C)" >>obtained; \
	  echo $(LCDD) $$file; \
	  ./$(LCDD) $$path | tr -d '\015' | sort | $(TEXT_MD5SUM) >>obtained; \
	done
	diff obtained $(srcdir)/$(EXPECTED)

check_cav18: $(LCDD)
	outfile=check_cav18_pplite; \
	echo $(TOUGH_EXAMPLES) >tough_examples; \
	echo $(ECHO_N) "$(ECHO_C)" >$$outfile; \
	export LC_ALL=C; \
	for path in $(srcdir)/examples/*.ext $(srcdir)/examples/*.ine; \
	do \
	  file=`basename $$path`; \
	  $(FGREP) $$file tough_examples >/dev/null 2>&1 && continue; \
	  echo $(ECHO_N) $$file " $(ECHO_C)" >>$$outfile; \
	  echo $(LCDD) $$path; \
	  ./$(LCDD) $$path | tr -d '\015' | sort | $(TEXT_MD5SUM) \
            >>$$outfile; \
	done

perf_cav18_all: $(LCDD)
	outfile=perf_cav18_all; \
	echo $(TOUGH_EXAMPLES) >tough_examples; \
	echo $(ECHO_N) "$(ECHO_C)" >$$outfile; \
	export LC_ALL=C; \
	for path in $(srcdir)/examples/*.ext $(srcdir)/examples/*.ine; \
	do \
	  file=`basename $$path`; \
	  $(FGREP) $$file tough_examples >/dev/null 2>&1 && continue; \
	  echo $(ECHO_N) $$file " $(ECHO_C)" >>$$outfile; \
	  echo $(LCDD) $$file "(3 times)"; \
	  timing=''; \
	  for i in 1 2 3; \
	  do \
	    singlerun=`./$(LCDD) -t $$path 2>&1 1>/dev/null | cut -d' ' -f 2`; \
	    echo $$singlerun; \
	    timing+=$$singlerun; \
	    timing+=' '; \
	  done; \
	  echo $$timing >>$$outfile; \
	done

CAV18_TAB3_EXAMPLES = \
cp6.ext \
cross12.ine \
in7.ine \
kkd38_6.ine \
kq20_11_m.ine \
metric80_16.ine \
mit31-20.ine \
mp6.ine \
reg600-5_m.ext \
sampleh8.ine \
trunc10.ine

perf_cav18_tab3: $(LCDD)
	outfile=perf_cav18_tab3; \
	echo $(CAV18_TAB3_EXAMPLES) >cav18_tab3_examples; \
	echo $(ECHO_N) "$(ECHO_C)" >$$outfile; \
	export LC_ALL=C; \
	for path in $(srcdir)/examples/*.ext $(srcdir)/examples/*.ine; \
	do \
	  file=`basename $$path`; \
	  $(FGREP) -v $$file cav18_tab3_examples >/dev/null 2>&1 && continue; \
	  echo $(ECHO_N) $$file " $(ECHO_C)" >>$$outfile; \
	  echo $(LCDD) $$file "(3 times)"; \
	  timing=''; \
	  for i in 1 2 3; \
	  do \
	    singlerun=`./$(LCDD) -t $$path 2>&1 1>/dev/null | cut -d' ' -f 2`; \
	    echo $$singlerun; \
	    timing+=$$singlerun; \
	    timing+=' '; \
	  done; \
	  echo $$timing >>$$outfile; \
	done


MOSTLYCLEANFILES = \
cav18_tab3_examples \
check_cav18_pplite \
obtained \
perf_cav18_all \
perf_cav18_tab3 \
tough_examples

endif HAVE_MD5SUM

$(top_builddir)/src/libpplite.la:
	$(MAKE) $(AM_MAKEFLAGS) -C $(top_builddir)/src libpplite.la
