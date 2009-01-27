dnl define(diversion_number, divnum)dnl
dnl divert(-1)


dnl Find the Fortran integer type that's equivalent to the C pointer type.
dnl
AC_DEFUN([UC_FORTRAN_PTR], [dnl
    AC_REQUIRE([UC_PROG_CC])
    AC_MSG_CHECKING(for Fortran integer type equivalent to C pointer)
    cat << EOF >conftest.c
#include <stdio.h>
int main()
{
    printf("%d\n", sizeof(char*));
    exit(0);
}
EOF
    doit='$CC -o conftest ${CFLAGS-} ${CPPFLAGS-} ${LDFLAGS-} conftest.c ${LIBS-}'
    if AC_TRY_EVAL(doit); then
        if type='integer*'`./conftest`; then
            AC_MSG_RESULT($type)
	    $1=$type
            AC_SUBST($1)
            unset type
        else
            AC_MSG_ERROR(Test program execution failure)
        fi
    else
        AC_MSG_ERROR(Test program build failure)
    fi
    rm conftest*
])


dnl Locate the perl(1) utility.
dnl
AC_DEFUN([UL_PROG_PERL], [dnl
    case "${PERL-unset}" in
        unset)
            AC_PROGRAMS_CHECK(PERL, perl)dnl
            case "$PERL" in
                '') UC_NEED_VALUE(PERL, [perl utility], /usr/local/bin/perl)
                    ;;
                *)  AC_SUBST(PERL)
                    ;;
            esac
            ;;
        *)
            AC_MSG_CHECKING(for perl utility)
            AC_SUBST(PERL)
            AC_MSG_RESULT($PERL)
            ;;
    esac
])


dnl Adjust compilation flags as necessary to achieve a position-independent
dnl extension library, if necessary.
AC_DEFUN([UL_PIC], [dnl
    AC_REQUIRE([UL_PERL_LINKTYPE])
    AC_MSG_CHECKING(for position-independent compilation flags)
    picflag=
    case `uname -s` in
        HP-UX)
            picflag=+z
            ;;
    esac
    case "$picflag" in
        '')
            AC_MSG_RESULT('')
            ;;
        *)
            AC_MSG_RESULT($picflag)
            UC_ENSURE(CFLAGS, $picflag)dnl
            UC_ENSURE(FFLAGS, $picflag)dnl
            UC_ENSURE(CXXFLAGS, $picflag)dnl
            ;;
    esac
])


dnl Determine the type of perl executable to create.
dnl
AC_DEFUN([UL_PERL_LINKTYPE], [dnl
    AC_MSG_CHECKING(for type of perl executable to create)
    case ${LINKTYPE-unset} in
        unset)
            case `uname -s` in
                ULTRIX)
                    LINKTYPE=static
                    ;;
                *)
                    LINKTYPE=dynamic
                    ;;
            esac
            ;;
        esac
    AC_MSG_RESULT($LINKTYPE)
    AC_SUBST(LINKTYPE)
])


dnl Determine the parameters for the top-level makefile regarding the
dnl perl/ subdirectory.
AC_DEFUN([UL_SUBDIR_PERL], [dnl
    AC_REQUIRE([UL_PROG_PERL])dnl
    AC_REQUIRE([UL_PERL_LINKTYPE])dnl
    case "${PERL}" in
        '')
            PERL_ALL=
            PERL_TEST=
            PERL_INSTALL=
            PERL_CLEAN=
            PERL_DISTCLEAN=
            PERL_MANUAL=
            ;;
        *)
            case "$LINKTYPE" in
                static)
                    PERL_ALL=perl/perl
                    PERL_INSTALL=perl/inst_perl
                    ;;
                *)
                    PERL_ALL=perl/dynamic
                    PERL_INSTALL=perl/install
                    ;;
            esac
            PERL_TEST=perl/test
            PERL_CLEAN=perl/clean
            PERL_DISTCLEAN=perl/distclean
            PERL_MANUAL=udunitsperl.1
            ;;
    esac
    AC_SUBST(PERL_ALL)dnl
    AC_SUBST(PERL_TEST)dnl
    AC_SUBST(PERL_INSTALL)dnl
    AC_SUBST(PERL_CLEAN)dnl
    AC_SUBST(PERL_DISTCLEAN)dnl
    AC_SUBST(PERL_MANUAL)dnl
])


dnl divert(diversion_number)dnl
