## $Rev:$ $Date:$
AC_INIT([GeneOut], [0.3], [dchaws+GeneOut@gmail.com])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_LANG(C++)
AC_PROG_CXX
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
 Makefile
 src/Makefile
])
AC_ARG_WITH(lapack-prefix,
[  --with-lapack-prefix=DIR   path to LAPACK],[
  if test "$withval" != "no"; then
    CPPFLAGS="${CPPFLAGS} -I$withval/include"
    LDFLAGS="${LDFLAGS} -L$withval/lib"
  fi
])

AC_ARG_WITH(blas-prefix,
[  --with-blas-prefix=DIR   path to BLAS],[
  if test "$withval" != "no"; then
    CPPFLAGS="${CPPFLAGS} -I$withval/include"
    LDFLAGS="${LDFLAGS} -L$withval/lib"
  fi
])

AC_ARG_WITH(svm-prefix,
[  --with-svm-prefix=DIR   path to SVM],[
  if test "$withval" != "no"; then
    CPPFLAGS="${CPPFLAGS} -I$withval/include"
    LDFLAGS="${LDFLAGS} -L$withval/lib"
  fi
])

AC_ARG_WITH([gmp],
  [AS_HELP_STRING([--with-gmp],
    [GMP Support. Optionally set to the location of gmp, or no to disable.])],
  [if test "$withval" != ""; then
    CPPFLAGS="${CPPFLAGS} -I$withval/include"
    LDFLAGS="${LDFLAGS} -L$withval/lib"
  fi
  ],
  [with_gmp=check])

AS_IF([test "x$with_gmp" != xno],
  [AC_CHECK_LIB([gmp], [__gmpz_init], ,
    [if test "x$with_gmp" != xcheck; then
       AC_MSG_FAILURE(
         [--with-gmp was given, but test for gmp failed])
     fi
    ])
   AC_CHECK_LIB([gmpxx], [main], ,
    [if test "x$with_gmp" != xcheck; then
       AC_MSG_FAILURE([--with-gmp was given, but test for gmpxx failed])
     fi
    ])
  ])



AC_CHECK_LIB(lapack, dgesvd_, ,
       [AC_MSG_ERROR([Lapack not found, see ])])
AC_CHECK_LIB(blas, main, ,
       [AC_MSG_ERROR([Blas not found, see ])])
AC_CHECK_LIB(svm, main, ,
       [AC_MSG_ERROR([SVM not found, see ])])
AC_OUTPUT

