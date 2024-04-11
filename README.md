# PPLite: convex polyhedra library for Abstract Interpretation

PPLite is an open-source C++ library
implementing the abstract domain of convex polyhedra,
to be used in tools for static analysis and verification.

<strong>Note:</strong> this public repository also collects
previous (non-versioned) releases of the library.

<h3>Current version
(see <a href="#available-downloads">below</a> for older ones)</h3>

2024-04-11:
<a href="releases/pplite-0.12.tar.gz">PPLite 0.12 can be downloaded</a>.

<h3>Support</h3>

If you need help for using PPLite or have some feature requests,
send an email to enea.zaffanella@unipr.it

<h3>Goals (and non-goals)</h3>

While being derived from the PPL (Parma Polyhedra Library), PPLite has
a very different goal: to provide researchers and students with a
lighter framework for experimenting with new ideas and algorithms in
the context of polyhedral computations.
In particular, PPLite is not aimed at implementing the full range of
abstract domains and operators made available by the PPL. The main
characteristics of PPLite are the following.
<ul>
  <li>Both closed and NNC rational convex polyhedra are supported.</li>
  <li>Performance and portability are deemed important, but not
      the main concern (ease of implementation and readability
      are given priority).</li>
  <li>The library is written in C++ (currently, c++17).</li>
  <li>The library is meant to be lightweight from the point of view
      of the developers: the goal is to reduce maintenance costs.
      This implies, among other things:
      <ul>
        <li>backward compatibility is a non-goal;
        <li>documentation is kept minimal (if not omitted altogether);
        <li>there is no plan to provide foreign language interfaces
            (but the library can be accessed through Apron,
             using the C, Java or OCaml bindings);
        <li>the library typically provides narrow contracts:
            preconditions on operators are not checked at runtime,
            so that user errors will lead to undefined behavior
            (rather than exceptions);
        <li>encapsulation is not fully enforced: a knowledgeable user
            can directly change the contents of inner data structures,
            e.g., to experiment with alternative implementations
            of domain operators.
      </ul>
</ul>


<h3>Compilation/installation instructions</h3>

Library development is based on Linux (Ubuntu).
The build system uses GNU autotools; the main external dependency
is the FLINT library, for arbitrary precision integer and rationals;
the library is written in C++ and requires a reasonably recent
C++ compiler (tested using g++ and clang++).

Installation of main requirements:

```
sudo apt-get install make autoconf automake libtool
sudo apt-get install libgmp-dev libmpfr-dev libflint-dev
```

Out-of-tree builds are strongly recommended:
in the following we assume using a source directory
named `MySources` and a build directory named `MyBuild`;
the commands shown should be adapted accordingly.

After cloning the library sources from git,
resulting in directory `MySources/PPLite`,
the following command will generate the configuration script
(note: this preliminary step is **not** needed if the sources
are extracted from a `.tar.gz` distribution):

```
cd MySources/PPLite
autoreconf --install
```

We then move to the build directory and configure it:

```
cd MyBuild
../MySources/PPLite/configure
```

This will use the default configuration options
(which can be overridden, see `configure --help`),
meaning that:
  * the installation prefix is `/usr/local`
  * the build is for release: optimizations are turned on
    and assertions are turned off
  * dependencies will be searched for in "standard" places

Use option `--prefix=INSTALLDIR` to change the installation directory.
Use option `--enable-assertions`
(and maybe also option `--enable-optimization=zero`) for a debug build.
If you have FLINT installed in a non-standard place the configure script
will complain: in that case you can use option
`--with-flint=DIR` or, if needed, the two options
`--with-flint-include=INCDIR` and `--with-flint-lib=LIBDIR`
to specify the location of the header files and the library;
similar options are available for MPFR and GMP libraries.

After configuration, the following commands will build the library,
run its tests and install it:

```
cd MyBuild
make -j8
make check -j8
sudo make install
```
Option `-j<N>` enables parallel builds using at most N workers;
consider removing the `sudo` prefix from the installation command
if it is not really needed.
Note that the (static and dynamic) library files will be installed
in `<prefix>/lib`, while the header files will be installed
in their own subdirectory `<prefix>/include/pplite`.


<h3>Main developers</h3>

The initial development team was composed by Anna Becchi and Enea Zaffanella.
Currently the library is maintained by Enea Zaffanella.

<h3>Contributors and past members of development team</h3>

Students collaborate to the development of the library (in a broad
sense, which is sometimes different from writing source code), under
supervision of Enea Zaffanella, during their internship and/or the
activities related to the final exam of the Bachelor and/or Master
Degree in Computer Science.

For the list of collaborators, see file CREDITS.

<A NAME="downloads">
<h3>Available downloads</h3>

<ul>

<li>
2024-04-11:
<a href="releases/pplite-0.12.tar.gz">PPLite 0.12 can be downloaded</a>.
<br>
The library no longer depends on the C++ interface of GMP (this change
mainly affects classes `pplite::Integer` and `pplite::Rational`).
The support for conditional thread safety is now enabled by default.
</li>

<li>
2023-06-10:
<a href="releases/pplite-0.11.tar.gz">PPLite 0.11 can be downloaded</a>.
<br>
While the library is based on the c++17 standard, its header files
are now c++11 standard compliant (hence can be used more easily in
other projects, e.g., the Apron library wrapper).
Added new method
  <code>Index_Set get_unconstrained() const;</code>
to polyhedra domains, returning the set of unconstrained space dimensions.
</li>

<li>
2023-05-29:
<strike><a href="releases/pplite-0.10.2.tar.gz">PPLite 0.10.2 can be downloaded</a>.</strike>
(<b>No longer maintained: switch to the most recent one.</b>)
<br>
Fixed a bug in the polyhedra simplification procedure.
The bug could be triggered when using NNC polyhedra and converting
from generator to constraints.
</li>

<li>
2023-03-07:
<strike><a href="releases/pplite-0.10.1.tar.gz">PPLite 0.10.1 can be downloaded</a>.</strike>
(<b>No longer maintained: switch to the most recent one.</b>)
<br>
Added support for the integral split operator.
</li>

<li>
2023-01-26:
<strike><a href="releases/pplite-0.10.tar.gz">PPLite 0.10 can be downloaded</a>.</strike>
(<b>No longer maintained: switch to the most recent one.</b>)
<br>
Main change: the Apron wrapper is no longer part of the library; it will
be integrated into Apron itself.
</li>

<li>
2023-01-02:
<strike><a href="releases/pplite-0.9.tar.gz">PPLite 0.9 can be downloaded</a>.</strike>
(<b>No longer maintained: switch to the most recent one.</b>)
<br>
The finite powerset domain is now a class template: pre-generated instances, also available via Apron interface, include finite sets of boxed polyhedra (P_Set) and finite sets of Cartesian factored boxed polyhedra (FP_Set).
</li>

<li>
2022-11-02:
<strike><a href="releases/pplite-0.8.tar.gz">PPLite 0.8 can be downloaded</a>.</strike>
(<b>No longer maintained: switch to the most recent one.</b>)
<br>
This new version adds a prototype implementation of the finite powerset
of Poly elements (PSet). It also provide an efficiency improved version
of the F_Poly domain and a couple of bug fixes in the Apron wrapper.
</li>

<li>
2022-08-26:
<strike><a href="releases/pplite-0.7.1.tar.gz">PPLite 0.7.1 can be downloaded</a>.</strike>
(<b>No longer maintained: don't use it, switch to the most recent one.</b>)
<br>
Note: this version fixes a silly build error; no functionality changes wrt 0.7.
</li>

<li>
2020-11-24:
<strike><a href="releases/pplite-0.7.tar.gz">PPLite 0.7 can be downloaded</a>.</strike>
(<b>No longer maintained: don't use it, switch to the most recent one.</b>)
<br>
New configuration option <em>--enable-apron</em> allows for compiling
the wrapper for (the C language bindings of) the Apron interface.
This version also adds a C++ polymorphic interface allowing to experiment
with several variants of the domain of convex polyhedra:
Poly, U_Poly, F_Poly, UF_Poly and their XXX_Stats versions,
computing timing information for abstract operators.
</li>

<li>
2020-04-23:
<strike><a href="releases/pplite-0.6.tar.gz">PPLite 0.6 can be downloaded</a>.</strike>
(<b>No longer maintained: don't use it, switch to the most recent one.</b>)
<br>
Added pplite::F_Poly, adding support for factored polyhedra
(experimental feature, not optimized).
</li>

<li>
2019-10-15:
<strike><a href="releases/pplite-0.5.1.tar.gz">PPLite 0.5.1 can be downloaded</a>.</strike>
(<b>No longer maintained: don't use it, switch to the most recent one.</b>)
<br>
This version fixes a bug affecting Poly::parallel_affine_image()
in version 0.5 (the bug was affecting computations of affine images
when the denominators of the expressions where different from 1).
</li>

<li>
2019-07-12:
<strike><a href="releases/pplite-0.5.tar.gz">PPLite 0.5 can be downloaded</a></strike>
(<b>No longer maintained: don't use it, switch to the most recent one.</b>)
<br>
Added helper class BBox to encode the (topologically closed,
possibly infinite) bounding box of a polyhedron;
this is used in new methods Poly::boxed_contains() and
Poly::boxed_is_disjoint_from() so as to speed up computation
(in particular, in client code computing over <em>sets</em> of polyhedra,
as found in tools such as PHAVerLite).
</li>

<li>
2019-03-08:
<strike><a href="releases/pplite-0.4.1.tar.gz">PPLite 0.4.1 can be downloaded</a></strike>
(<b>No longer maintained: don't use it, switch to the most recent one.</b>)
<br>
(Note: version 0.4.1 fixes a bug affecting Poly::hash() in version 0.4,
which was released on 2019-03-03.)
Added new method Poly::hash(), which can be used to quickly decide
that two polyhedra are different. Added new method Poly::con_hull_assign(),
computing the <em>constraint hull</em> of two polyhedra (note: the result
depends on the specific constraint representation computed by the library).
Improved the efficiency of the parallel affine image method.
Corrected a few, corner case bugs (affected functionalities were:
computation of the bhrz03 widening; the removal of space dimensions;
the addition of generating lines).
</li>

<li>
2018-10-09:
<strike><a href="releases/pplite-0.3.tar.gz">PPLite 0.3 can be downloaded</a>.</strike>
(<b>No longer maintained: don't use it, switch to the most recent one.</b>)
<br>
Added wrappers on the polyhedra domain:
Poly_Stats collects time statistics on abstract operator calls;
U_Poly improves the handling of unconstrained space dimensions.
Corrected a few bugs: intersection of 0-dim polyhedra;
printing of universe polyhedra; minimization of an affine function;
a couple of other bugs in code used only in debugging mode.
</li>

<li>
2018-07-10:
<strike><a href="releases/pplite-0.2.tar.gz">PPLite 0.2 can be downloaded</a>.</strike>
(<b>No longer maintained: don't use it, switch to the most recent one.</b>)
<br>
Added support for (conditional) thread safety:
it is now possible to develop multithreaded applications using PPLite.
Implemented a few other operators on the polyhedra domain,
including the more precise bhrz03 widening.
</li>

<li>
2018-05-08:
<strike><a href="releases/pplite-0.1.tar.gz">PPLite 0.1 can be downloaded</a>.</strike>
(<b>No longer maintained: don't use it, switch to the most recent one.</b>)
<br>
Initial release.
</li>

</ul>
