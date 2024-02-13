AC_INIT([PPLite],[0.12],[enea.zaffanella@unipr.it],[pplite],[github.com/ezaffanella/PPLite])

# Specify the location of additional local Autoconf macros.
AC_CONFIG_MACRO_DIR([m4])

AC_CANONICAL_HOST

AM_INIT_AUTOMAKE([foreign no-define silent-rules 1.11])

# Generate a configuration header file.
AC_CONFIG_HEADERS([config.h])

AH_BOTTOM([
#if defined(PPLITE_NDEBUG) && !defined(NDEBUG)
# define NDEBUG PPLITE_NDEBUG
#endif
])

# Rename configuration header file (and prefix its defines).
AX_PREFIX_CONFIG_H(pplite-config.h)

ISODATE=`date +%Y-%m-%d`
AC_SUBST(ISODATE)

# C compiler.

AC_ARG_WITH(cc,
  AS_HELP_STRING([--with-cc=XXX], [use XXX as the C compiler]),
  CC=$with_cc)

# C++ compiler.

AC_ARG_WITH(cxx,
  AS_HELP_STRING([--with-cxx=XXX], [use XXX as the C++ compiler]),
  CXX=$with_cxx)

# Checks for programs.

save_CFLAGS="$CFLAGS"
AC_PROG_CC
CFLAGS="$save_CFLAGS"

AM_CONDITIONAL(GCC, test x"$GCC" = xyes)

save_CXXFLAGS="$CXXFLAGS"
AC_PROG_CXX
CXXFLAGS="$save_CXXFLAGS"

# Check if the C++ compiler supports C++17 features.
AX_CXX_COMPILE_STDCXX(17, noext, mandatory)

AC_PROG_CXXCPP

# Compilation flags.

COMP_FLAGS=""
OPT_FLAGS=""

# --enable-debugging
enableval=yes
AC_MSG_CHECKING([whether to compile with debug info])
AC_ARG_ENABLE(debugging,
  AS_HELP_STRING([--enable-debugging], [compile with debugging information]))
case "${enableval}" in
yes)
  AC_MSG_RESULT(yes)
  enable_debugging=yes
  OPT_FLAGS="-g"
  ;;
no)
  AC_MSG_RESULT(no)
  enable_debugging=no
  ;;
*)
  AC_MSG_ERROR([bad value ${enableval} for --enable-debugging, needs yes or no])
  enable_debugging=yes
  ;;
esac
AM_CONDITIONAL(DEBUGGING_ENABLED, test x"$enable_debugging" = xyes)

# --enable-profiling
enableval=no
AC_MSG_CHECKING([whether to compile for profiling])
AC_ARG_ENABLE(profiling,
  AS_HELP_STRING([--enable-profiling], [compile for profiling]))
case "${enableval}" in
yes)
  AC_MSG_RESULT(yes)
  OPT_FLAGS="-g"
  COMP_FLAGS="$COMP_FLAGS -pg"
  ;;
no)
  AC_MSG_RESULT(no)
  ;;
*)
  AC_MSG_ERROR([bad value ${enableval} for --enable-profiling, needs yes or no])
  ;;
esac

# --enable-coverage
enableval=no
AC_MSG_CHECKING([whether to compile for test coverage])
AC_ARG_ENABLE(coverage,
  AS_HELP_STRING([--enable-coverage], [compile for test coverage]))
case "${enableval}" in
yes)
  AC_MSG_RESULT(yes)
  OPT_FLAGS="-g"
  COMP_FLAGS="$COMP_FLAGS -fprofile-arcs -ftest-coverage"
  extra_libraries="${extra_libraries} -lgcov"
  ;;
no)
  AC_MSG_RESULT(no)
  ;;
*)
  AC_MSG_ERROR([bad value ${enableval} for --enable-coverage, needs yes or no])
  ;;
esac

# --enable-assertions
enableval=no
AC_MSG_CHECKING([whether to enable checking of run-time assertions])
AC_ARG_ENABLE(assertions,
  AS_HELP_STRING([--enable-assertions], [check run-time assertions]))
case "${enableval}" in
yes)
  AC_MSG_RESULT(yes)
  enable_assertions=yes
  ;;
no)
  AC_MSG_RESULT(no)
  enable_assertions=no
  AC_DEFINE(NDEBUG, 1, [Assertions are disabled when this is defined.])
  debug_flag="-DNDEBUG=1"
  ;;
*)
  AC_MSG_ERROR([bad value ${enableval} for --enable-assertions, needs yes or no])
  enable_assertions=no
  ;;
esac
AM_CONDITIONAL(ASSERTIONS_ENABLED, test x"$enable_assertions" = xyes)

# --enable-valgrind-tests
enableval=no
AC_MSG_CHECKING([whether tests should be run under Valgrind])
AC_ARG_ENABLE(valgrind-tests,
  AS_HELP_STRING([--enable-valgrind-tests],
                 [run library tests under Valgrind]))
case "${enableval}" in
yes)
  AC_CHECK_PROG(VALGRIND, valgrind, valgrind)
  if test x"$VALGRIND" = xvalgrind
  then
    AC_MSG_RESULT(yes)
    valgrind_tests=yes
  else
    AC_MSG_ERROR([bad value ${enableval} for --enable-valgrind-tests, cannot find Valgrind])
  fi
  ;;
no)
  AC_MSG_RESULT(no)
  valgrind_tests=no
  ;;
*)
  AC_MSG_ERROR([bad value ${enableval} for --enable-valgrind-tests, needs yes or no])
  ;;
esac
AM_CONDITIONAL(VALGRIND_TESTS_ENABLED, test x"$valgrind_tests" = xyes)

# --enable-optimization
arch=no
enableval=standard
AC_MSG_CHECKING([whether to enable optimizations])
AC_ARG_ENABLE(optimization,
  AS_HELP_STRING([--enable-optimization@<:@=LEVEL@:>@],
                 [enable compiler optimizations]))
case "${enableval}" in
sspeed)
  AC_MSG_RESULT(sspeed)
  OPT_FLAGS="$OPT_FLAGS -O3 -fomit-frame-pointer"
  arch=yes
  ;;
speed)
  AC_MSG_RESULT(speed)
  OPT_FLAGS="$OPT_FLAGS -O3"
  arch=yes
  ;;
size)
  AC_MSG_RESULT(size)
  OPT_FLAGS="$OPT_FLAGS -Os"
  arch=yes
  ;;
standard | yes)
  AC_MSG_RESULT(standard)
  OPT_FLAGS="$OPT_FLAGS -O2"
  ;;
mild)
  AC_MSG_RESULT(mild)
  OPT_FLAGS="$OPT_FLAGS -O1"
  ;;
zero)
  AC_MSG_RESULT(zero)
  OPT_FLAGS="$OPT_FLAGS -O0"
  ;;
no)
  AC_MSG_RESULT(no)
  OPT_FLAGS="$OPT_FLAGS -O0"
  ;;
*)
  AC_MSG_ERROR([bad value ${enableval} for --enable-optimization, needs sspeed, speed, size, standard, mild, zero, yes or no])
  ;;
esac

# --enable-arch
enableval=$arch
AC_MSG_CHECKING([for which architecture to optimize])
AC_ARG_ENABLE(arch,
  AS_HELP_STRING([--enable-arch@<:@=ARCH@:>@],
                 [optimize for architecture ARCH]))
case "${enableval}" in
yes)
  m=`uname -m`
  case $m in
  i?86 | k6 | athlon)
    AC_MSG_RESULT($m)
    OPT_FLAGS="$OPT_FLAGS -march=$m"
    ;;
  *)
    AC_MSG_RESULT(default)
    ;;
  esac
  ;;
no)
  AC_MSG_RESULT(default)
  ;;
*)
  AC_MSG_RESULT($enableval)
  OPT_FLAGS="$OPT_FLAGS -march=$enableval"
  ;;
esac

# --enable-thread-safe
enableval=no
AC_MSG_CHECKING([whether to compile for thread safety])
AC_ARG_ENABLE(thread-safe,
  AS_HELP_STRING([--enable-thread-safe], [compile for thread safety]))
case "${enableval}" in
yes)
  AC_MSG_RESULT(yes)
  enable_thread_safe=yes
  AC_DEFINE(THREAD_SAFE, 1, [The library is (conditionally) thread safe when this is defined.])
  COMP_FLAGS="$COMP_FLAGS -pthread"
  extra_libraries="${extra_libraries} -pthread"
  ;;
no)
  AC_MSG_RESULT(no)
  enable_thread_safe=no
  ;;
*)
  AC_MSG_ERROR([bad value ${enableval} for --enable-thread-safe, needs yes or no])
  ;;
esac
AM_CONDITIONAL(THREAD_SAFE, test x"$enable_thread_safe" = xyes)


CFLAGS="$COMP_FLAGS $OPT_FLAGS $CFLAGS"
CXXFLAGS="$COMP_FLAGS $OPT_FLAGS $CXXFLAGS"

# Allow additions to C compilation flags.
AC_ARG_WITH(cflags,
  AS_HELP_STRING([--with-cflags=XXX],
                 [add XXX to the options for the C compiler]),
  CFLAGS="$CFLAGS $with_cflags")

# Allow additions to C++ compilation flags.
AC_ARG_WITH(cxxflags,
  AS_HELP_STRING([--with-cxxflags=XXX],
                 [add XXX to the options for the C++ compiler]),
  CXXFLAGS="$CXXFLAGS $with_cxxflags")

# Libtool.
LT_PREREQ([1.5.24])
LT_INIT([dlopen])
AM_CONDITIONAL(ENABLE_SHARED, test x"$enable_shared" = xyes)
AM_CONDITIONAL(ENABLE_STATIC, test x"$enable_static" = xyes)

# Make sure libtool is used for linking.
AC_USE_LIBTOOL

# --enable-integers
# Checks whether to use Flint (default) or GMP integers.
enableval=flint
AC_MSG_CHECKING([the type of unbounded integers to use as coefficients])
AC_ARG_ENABLE(integers,
  AS_HELP_STRING([--enable-integers=TYPE],
                 [select the type of the integers (flint or gmp)]))
case "${enableval}" in
flint)
  AC_MSG_RESULT([FLINT fmpz])
  integer_kind=flint
  ;;
gmp)
  AC_MSG_RESULT([GMP mpz])
  integer_kind=gmp
  ;;
*)
  AC_MSG_ERROR([bad value ${enableval} for --enable-integers, flint or gmp])
  ;;
esac
AM_CONDITIONAL(USE_FLINT_INTEGERS, test x"$integer_kind" = xflint)

if test x"$integer_kind" = xflint
then
  AC_DEFINE(USE_FLINT_INTEGERS, 1, [Defined if the integral type used to represent coefficients is FLINT's fmpz.])
else
  AC_DEFINE(USE_FLINT_INTEGERS, 0, [Defined if the integral type used to represent coefficients is FLINT's fmpz.])
fi

# Checks on the GMP library.
AC_CHECK_GMP
if test x"$have_gmp" = xno
then
  AC_MSG_ERROR([Cannot find GMP version 4.1.3 or higher.])
fi

extra_includes="${extra_includes}${extra_includes:+ }${gmp_include_options}"
extra_libraries="${extra_libraries}${extra_libraries:+ }${gmp_library_options}"

if test x"$integer_kind" = xflint
then
  # Checks on the MPFR library.
  AC_CHECK_MPFR
  # Checks on the Flint library.
  AC_CHECK_FLINT
  if test x"$have_mpfr" = xno
  then
    AC_MSG_ERROR([Cannot find MPFR.])
  elif test x"$have_flint" = xno
  then
    AC_MSG_ERROR([Cannot find Flint version 2.5.2 or higher.])
  else
    extra_includes="${extra_includes}${extra_includes:+ }${flint_include_options} ${mpfr_include_options}"
    extra_libraries="${extra_libraries}${extra_libraries:+ }${flint_library_options} ${mpfr_library_options}"
  fi
fi

# Checks for md5sum.
AC_CHECK_PROG(MD5SUM, md5sum, md5sum)
AM_CONDITIONAL(HAVE_MD5SUM, test x$MD5SUM = xmd5sum)
if test x$MD5SUM = xmd5sum
then
  AC_TEXT_MD5SUM
fi

# We want to compile with lots of warnings enabled.
if test x"$GCC" = xyes
then
  CFLAGS="$CFLAGS -Wall -Wextra"
fi
if test x"$GXX" = xyes
then
  CXXFLAGS="$CXXFLAGS -Wall -Wextra"
fi

AC_SUBST(debug_flag)
AC_SUBST(extra_includes)
AC_SUBST(extra_libraries)

AC_CONFIG_FILES(Makefile
                src/Makefile
                tests/Makefile
                tests/Con/Makefile
                tests/Gen/Makefile
                tests/Index_Set/Makefile
                tests/Integer/Makefile
                tests/Linear_Expr/Makefile
                tests/F_Poly/Makefile
                tests/Poly/Makefile
                tests/Poly_Impl/Makefile
                tests/PSet/Makefile
                tests/Var/Makefile
                demos/Makefile
                demos/pplite_lcdd/Makefile
                demos/pplite_lcdd/examples/Makefile
                )
AC_OUTPUT
