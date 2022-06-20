# PPLite: convex polyhedra library for Abstract Interpretation

PPLite is an open-source C++ library implementing the abstract domain of convex polyhedra, to be used in tools for static analysis and verification.

<h3>Current version (see below for older ones)</h3>

    2020-11-24: PPLite 0.7 can be downloaded.

<h3>Support</h3>

  If you need help for using PPLite or have some feature requests, send an email to enea.zaffanella@unipr.it

<h3>Goals (and non-goals)</h3>

While being derived from the PPL (Parma Polyhedra Library), PPLite has a very different goal: to provide researchers and students with a lighter framework for experimenting with new ideas and algorithms in the context of polyhedral computations. In particular, PPLite is not aimed at implementing the full range of abstract domains and operators made available by the PPL. The main characteristics of PPLite are the following.
<ul>
  <li>Both closed and NNC rational convex polyhedra are supported.</li>
  <li>Exact computations are based on FLINT</li>
  <li>Performance is deemed important, but not the main concern
     (ease of implementation and readability are given priority).</li>
  <li>Portability is deemed important, but not the main concern
      (ease of implementation and readability are given priority).</li>
  <li>The library is written in modern C++. The developers should feel free to use language features
      made available by the recent standards, provided these lead to a simpler implementation
    or improve on code readability.</li>
  <li>The library is meant to be lightweight from the point of view of the developers:
      the goal is to reduce maintenance costs. This implies, among other things:
      <ul>
        <li>backward compatibility is a non-goal;
        <li>documentation is kept minimal (if not omitted altogether);
        <li>there is no plan to provide foreign language interfaces;
        <li>the library typically provides narrow contracts: preconditions on operators are not checked at runtime, so that user errors will lead to undefined behavior (rather than exceptions);
        <li>encapsulation is not fully enforced: a knowledgeable user can directly change the contents of inner data structures, e.g., to experiment with alternative implementations of domain operators.
      </ul>
</ul>

<h3>Main developers</h3>

    Enea Zaffanella (supervisor)
    Anna Becchi (former student)

<h3>Contributors and past members of development team</h3>

    Gino Ceresini (former student)
    Carlotta Colla (former student)
    Maria Chiara Colla (former student)
    Fabio Cristini (former student)
    Eduard Ispas (former student)
    Lorenzo Mora (student)
    Riccardo Mori (student)
    Sara Musiari (student)
    Velia Pierdomenico (student)
    Dana Greta Pop (student)
    Luigi Zaccone (former student)


<A NAME="downloads">
<h3>Available downloads</h3>
<ul>
<li>
2020-11-24:
<a href="releases/pplite-0.7.tar.gz">PPLite 0.7 can be downloaded</a>.
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
(<b>No longer maintained: don't use it, switch to the most recent one</b>).
<br>
Added pplite::F_Poly, adding support for factored polyhedra
(experimental feature, not optimized).
</li>
<li>
2019-10-15:
<strike><a href="releases/pplite-0.5.1.tar.gz">PPLite 0.5.1 can be downloaded</a>.</strike>
(<b>No longer maintained: don't use it, switch to the most recent one</b>).
<br>
This version fixes a bug affecting Poly::parallel_affine_image()
in version 0.5 (the bug was affecting computations of affine images
when the denominators of the expressions where different from 1).
</li>
<li>
2019-07-12:
<strike><a href="releases/pplite-0.5.tar.gz">PPLite 0.5 can be downloaded</a></strike>
(<b>No longer maintained: don't use it, switch to the most recent one</b>).
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
<strike><a href="releases/pplite-0.4.1.tar.gz">PPLite 0.4.1 can be downloaded</a></strike> (<b>No longer maintained: don't use it, switch to the most recent one</b>).
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
<strike><a href="releases/pplite-0.3.tar.gz">PPLite 0.3 can be downloaded</a>.</strike> (<b>No longer maintained: don't use it, switch to the most recent one</b>.)
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
<strike><a href="releases/pplite-0.2.tar.gz">PPLite 0.2 can be downloaded</a>.</strike> (<b>No longer maintained: don't use it, switch to the most recent one</b>.)
<br>
Added support for (conditional) thread safety:
it is now possible to develop multithreaded applications using PPLite.
Implemented a few other operators on the polyhedra domain,
including the more precise bhrz03 widening.
</li>
<li>
2018-05-08:
<strike><a href="releases/pplite-0.1.tar.gz">PPLite 0.1 can be downloaded</a>.</strike> (<b>No longer maintained: don't use it, switch to the most recent one</b>.)
<br>
Initial release.
</li>
</ul>
