dnl Setup for making a manual-page database.
dnl
AC_DEFUN([UC_MAKEWHATIS],
[
    #
    # NB: We always want to define WHATIS to prevent the
    # $(MANDIR)/$(WHATIS) make(1) target from being just $(MANDIR)/ and
    # conflicting with the (directory creation) target with the same name.
    #
    WHATIS=whatis
    case `uname -sr` in
        (BSD/OS*)
            # Can't generate a user-database -- only /usr/share/man/whatis.db.
            MAKEWHATIS_CMD=
            ;;
        ('IRIX64 6.5'|'IRIX 6.5')
            MAKEWHATIS_CMD='/usr/lib/makewhatis -M $(MANDIR) $(MANDIR)/whatis'
            ;;
        ('IRIX 6'*)
            # Can't generate a user-database.
            MAKEWHATIS_CMD=
            ;;
        (HP-UX*)
            # Can't generate a user-database -- only /usr/lib/whatis.
            MAKEWHATIS_CMD=
            ;;
        ('Linux '*)
            # /usr/sbin/makewhatis doesn't work
            MAKEWHATIS_CMD=
            ;;
        (ULTRIX*)
            # Can't generate a user-database -- only /usr/lib/whatis.
            MAKEWHATIS_CMD=
            ;;
        (*)
            if test -r /usr/man/windex; then
                WHATIS=windex
            fi
            AC_CHECK_PROGS(prog, catman makewhatis /usr/lib/makewhatis)
            case "$prog" in
                (*catman*)
                    MAKEWHATIS_CMD=$prog' -w -M $(MANDIR)'
                    ;;
                (*makewhatis*)
                    MAKEWHATIS_CMD=$prog' $(MANDIR)'
                    ;;
            esac
            ;;
    esac
    AC_SUBST(WHATIS)
    AC_SUBST(MAKEWHATIS_CMD)
    AC_MSG_CHECKING(for manual-page index command)
    AC_MSG_RESULT($MAKEWHATIS_CMD)
])


dnl Initialize this Unidata autoconf(1)-support module without support
dnl for using the port/ subdirectory other than for the makefiles.
dnl
AC_DEFUN([UC_INIT_NOPORT],
[dnl
dnl
dnl Don't ask me why, but the ac_require() around uc_customize() is 
dnl necessary to prevent the uc_port stuff from being written to the
dnl output file before the uc_customize stuff.
dnl
    AC_REQUIRE([UC_CUSTOMIZE])
    rm -f confdefs.missing
    UC_MAKEFILE(port/master.mk) dnl
    UC_MAKEFILE(port/Makefile) dnl
    UC_DEFAULT(PORT_SUBDIRS, ) dnl
    UC_DEFAULT(PORT_MANIFEST, ) dnl
    AC_SUBST(PORT_ALL) dnl
    AC_SUBST(PORT_CLEAN) dnl
    AC_SUBST(PORT_DISTCLEAN) dnl
    AC_SUBST(PORT_INSTALL) dnl
    UC_MAKEWHATIS
])


dnl Initialize this Unidata autoconf(1)-support module.
dnl
AC_DEFUN([UC_INIT],
[dnl
    AC_REQUIRE([UC_INIT_NOPORT])
    AC_CONFIG_HEADER(port/misc/udposix.h)
    UC_PORT
])


dnl Define a substitution for both `@...@' and C preprocessor forms.
dnl
AC_DEFUN([UC_REPLACE], [dnl
$1="$2"
AC_SUBST($1)dnl
AC_DEFINE($1,$2)
])


dnl Set the value of a variable.  Use the environment if possible; otherwise
dnl set it to a default value.  Call the substitute routine.
dnl
AC_DEFUN([UC_DEFAULT], [dnl
$1=${$1-"$2"}
AC_SUBST([$1])
])


dnl Configure a C header file.  This differs from AC_CONFIG_HEADER by
dnl accumulating the specified headers into the list rather than redefining
dnl the list each time.
dnl
dnl UC_CONFIG_HEADER(output_header_file)
dnl
AC_DEFUN([UC_CONFIG_HEADER], [dnl
ifdef([AC_LIST_HEADER], , [define([AC_LIST_HEADER], [])]) dnl
ifelse(index(AC_LIST_HEADER, $1), -1, [dnl
dnl    echo Adding $1 to configuration headers
    define([tmp], AC_LIST_HEADER) dnl
    define([AC_LIST_HEADER], tmp $1) dnl
    undefine([tmp])] dnl
) dnl
dnl echo configuration headers: AC_LIST_HEADER
])


dnl Set up for customizing the port/ subdirectory.
dnl
AC_DEFUN([UC_PORT],
[dnl
    UC_UDPOSIX
    UC_PROG_AR
    UC_PROG_TAR
    AC_PROG_RANLIB
    UC_BINFTP
    UC_MAKEFILE(port/Makefile) dnl
    UC_MAKEFILE(port/master.mk) dnl
    UC_MAKEFILE(port/misc/Makefile) dnl
    UC_DEFAULT(UC_LIBOBJS, ) dnl
    UC_DEFAULT(PORT_HEADERS, ) dnl
    UC_DEFAULT(PORT_MANIFEST, ) dnl
    UC_DEFAULT(PORT_SUBDIRS, ) dnl

    UC_ENSURE(PORT_ALL, misc/all)
    UC_ENSURE(PORT_INSTALL, misc/install)
    UC_ENSURE(PORT_CLEAN, misc/clean)
    UC_ENSURE(PORT_DISTCLEAN, misc/distclean)

    AC_SUBST(PORT_ALL) dnl
    AC_SUBST(PORT_INSTALL) dnl
    AC_SUBST(PORT_CLEAN) dnl
    AC_SUBST(PORT_DISTCLEAN) dnl

    UC_DEFAULT(PORT_CPP_LDM, ) dnl
    udportdir=`pwd`/port/misc
    UC_LINK_REF(LD_UDPORT, $udportdir, udport) dnl
    AC_SUBST(LD_UDPORT)
    unset udportdir
])


dnl Finish with everything (both GNU and Unidata autoconf(1) support).
dnl
AC_DEFUN([UC_FINISH], 
[dnl
    UC_CHECK_MISSING dnl
    AC_OUTPUT($1 ${CONFIG_FILES-},
            [dnl
                UC_POSTPROCESS_MAKEFILES(${CONFIG_FILES-}) dnl
            ],
            [dnl
                CC="$CC"
            ]dnl
    ) dnl
])


dnl Add a member to the list of makefiles.
dnl
AC_DEFUN([UC_MAKEFILE], [dnl
UC_ENSURE(CONFIG_FILES, $1)dnl
])


dnl Determine the operating system.
dnl
AC_DEFUN([UC_OS], [dnl
    AC_MSG_CHECKING(type of operating system) 
    if test -z "$OPSYS"; then
      OPSYS=`uname -s | tr '[[A-Z]]' '[[a-z]]' | sed 's;/;;g'`
      if test -z "$OPSYS"; then
        UC_NEED_VALUE(OPSYS, [operating system], sunos5)dnl
      fi
    fi
    case $OPSYS in
        aix)
            OS_NAME=`uname -s`
            OS_MAJOR=`uname -v | sed 's/[[^0-9]]*\([[0-9]]*\)\..*/\1/'`
            ;;
        hp-ux)
            OPSYS=hpux`uname -r | sed 's/[[A-Z.0]]*\([[0-9]]*\).*/\1/'`
            OS_NAME=HPUX
            OS_MAJOR=`uname -r | sed 's/[[A-Z.0]]*\([[0-9]]*\).*/\1/'`
            ;;
        irix)
            OPSYS=${OPSYS}`uname -r | sed 's/\..*//'`
            OS_NAME=IRIX
            OS_MAJOR=`uname -r | sed 's/\..*//'`
            ;;
        osf*)
            OS_NAME=OSF1
            OS_MAJOR=`uname -r | sed 's/[[^0-9]]*\([[0-9]]*\)\..*/\1/'`
            ;;
        sn*)
            OPSYS=unicos
            OS_NAME=UNICOS
            OS_MAJOR=`uname -r | sed 's/[[^0-9]]*\([[0-9]]*\)\..*/\1/'`
            ;;
        sunos)
            OS_NAME=SunOS
            OS_MAJOR=`uname -r | sed 's/\..*//'`
            OPSYS=$OPSYS$OS_MAJOR
            ;;
        ultrix)
            OS_NAME=ULTRIX
            OS_MAJOR=`uname -r | sed 's/\..*//'`
            ;;
        *)
            # On at least one UNICOS system, 'uname -s' returned the
            # hostname (sigh).
            if uname -a | grep CRAY >/dev/null; then
                OPSYS=unicos
                OS_NAME=UNICOS
            else
                OS_NAME=`uname -s | sed 's/[[^A-Za-z0-9_]]//g'`
            fi
            OS_MAJOR=`uname -r | sed 's/[[^0-9]]*\([[0-9]]*\)\..*/\1/'`
            ;;
    esac

    # Adjust OPSYS for CRAY MPP environment.
    #
    case "$OPSYS" in
    unicos)
        AC_REQUIRE([UC_PROG_CC])
        case "$CC$TARGET$CFLAGS" in
        *cray-t3*)
            OPSYS=unicos-mpp
            ;;
        esac
        ;;
    esac

    AC_SUBST(OPSYS)dnl
    AC_DEFINE_UNQUOTED(OS_NAME, $OS_NAME)
    AC_DEFINE_UNQUOTED(OS_MAJOR, $OS_MAJOR)
    AC_MSG_RESULT($OPSYS) 
])


dnl Check for C++ compiler.  This macro is used rather than using the
dnl autoconf macro directly because we want it to obtain the native C++
dnl compiler if possible.
dnl
AC_DEFUN([UC_PROG_CXX], [dnl
    case ${CXX+set} in
        set)
            AC_MSG_CHECKING(for C++ compiler)
            AC_MSG_RESULT($CXX)
            AC_SUBST(CXX)dnl
            ;;
        *)
            case `uname -s` in
                AIX*)
                    native_cxx=xlC
                    ;;
                HP-UX*)
                    native_cxx=CC
                    ;;
                IRIX*)
                    native_cxx=CC
                    ;;
                OSF*)
                    native_cxx=cxx
                    ;;
                SunOS)
                    native_cxx=CC
                    ;;
                *)
                    native_cxx=
                    ;;
            esac
            AC_CHECK_PROGS(CXX, $native_cxx CC cxx c++ g++ gcc, )dnl
            case "${CXX-unset}" in
                unset)
                    UC_NEED_VALUE(CXX, [C++ compiler], CC)dnl
                    ;;
            esac
            ;;
    esac
])


dnl Check for C compiler.  This macro replaces the ac_prog_cc macro because 
dnl that macro prefers the GNU C compiler.
dnl
AC_DEFUN([UC_PROG_CC],
[dnl
dnl
dnl The AC_REQUIRE macro runs a diversion stack.  If it gets too deep, then
dnl some things are output in the wrong order.  Thus, we use AC_BEFORE here
dnl rather than AC_REQUIRE where we'd like.
    AC_BEFORE([$0], [UC_FINISH])
    AC_BEFORE([$0], [UC_CONST])
    AC_BEFORE([$0], [UC_VOLATILE])
    AC_BEFORE([$0], [UC_SIGNED])
    AC_BEFORE([$0], [UC_PROTOTYPES])
    AC_BEFORE([$0], [UC_MAKESTRING])
    AC_BEFORE([$0], [UC_GLUE])
    AC_BEFORE([$0], [UC_VOIDP])
    AC_BEFORE([$0], [UC_FUNC_DECL])
    AC_BEFORE([$0], [UC_TYPEDEF])
    AC_BEFORE([$0], [UC_STRUCT])
    AC_BEFORE([$0], [UC_MACRO])
    AC_BEFORE([$0], [UC_UTHREADS])
    AC_BEFORE([$0], [UC_UDPOSIX_FLOAT])
    AC_BEFORE([$0], [UC_PROG_CC_MAKEDEPEND])
    #
    # Ensure that the CC variable is unset so that it can be
    # set here rather than by the autoconf-generated
    # configure-script preamble.
    #
    # unset CC
    #
    case ${CC-unset} in
        unset)
            case `uname -s` in
                AIX)
                    AC_CHECK_PROGS(CC, c89 xlc cc gcc) dnl
                    ;;
                HP-UX)
                    AC_CHECK_PROGS(CC, c89 cc gcc) dnl
                    ;;
                IRIX*)
                    AC_CHECK_PROGS(CC, cc gcc) dnl
                    ;;
                OSF1|ULTRIX)
                    AC_CHECK_PROGS(CC, cc gcc) dnl
                    case "$CC" in
                    cc)
                        case `uname -m` in
                        VAX)
                            ;;
                        *)
                            UC_ENSURE(CPPFLAGS, -std) dnl
                            ;;
                        esac
                        ;;
                    esac
                    ;;
                SunOS)
                    case `uname -r` in
                        4*)
                            AC_CHECK_PROGS(CC, acc cc gcc) dnl
                            ;;
                        5*)
                            AC_CHECK_PROGS(CC, cc gcc) dnl
#
#                           The following is commented-out because
#                           the configure script uses CPPFLAGS when
#                           compiling C++ source and SunOS 5's CC (at
#                           least) emits error messages when given the
#                           -Xa option causing the configure script to
#                           abandon `$CXX -E' in favor of `/lib/cpp'.
#
#                           case "$CC" in
#                               *gcc*)
#                                   ;;
#                               *)
#                                   UC_ENSURE(CPPFLAGS, -Xa) dnl
#                                   ;;
#                           esac
                            ;;
                    esac
                    ;;
                *)
                    AC_CHECK_PROGS(CC, c89 cc gcc) dnl
                    ;;
            esac
            ;;
        *)
            AC_MSG_CHECKING(for C compiler)
            AC_MSG_RESULT($CC)
            ;;
    esac
    case "${CC-}" in
    '')
        UC_NEED_VALUE(CC, [C compiler], /bin/cc) dnl
        ;;
    *)
        # Find out if we are using GNU C, under whatever name.
        cat <<UD_EOF > conftest.c
#ifdef __GNUC__
            yes
#endif
UD_EOF
        ${CC} -E conftest.c > conftest.out 2>&1
        if egrep yes conftest.out >/dev/null 2>&1; then
            GCC=1 # For later tests.
        fi
        AC_SUBST(CC) dnl
        case `uname -s` in
            AIX)
                UC_ENSURE(CPPFLAGS, -D_ALL_SOURCE) dnl
                AC_DEFINE(_ALL_SOURCE) dnl
                ;;
            HP-UX)
                UC_ENSURE(CPPFLAGS, -D_HPUX_SOURCE) dnl
                AC_DEFINE(_HPUX_SOURCE) dnl
                ;;
        esac
        ;;
    esac
    rm -f conftest*
])


dnl Check for functioning `const' keyword
dnl
AC_DEFUN([UC_CONST],
[dnl
    AC_MSG_CHECKING(for C const)
    AC_TRY_COMPILE(
        ,
        [/* Ultrix mips cc rejects this.  */
            typedef int charset[2]; const charset x;
        ],
        [AC_MSG_RESULT(yes)],
        [
           UC_REPLACE(UD_NO_CONST,1)dnl
           AC_MSG_RESULT(no)
        ]
    )
])


dnl Check for functioning `volatile' keyword
dnl
AC_DEFUN([UC_VOLATILE], [dnl
    AC_MSG_CHECKING(for C volatile)
    AC_TRY_COMPILE(
        ,
        [typedef int charset[2]; volatile charset x;],
        [AC_MSG_RESULT(yes)],
        [
            AC_MSG_RESULT(no)
            UC_REPLACE(UD_NO_VOLATILE,1)
        ]
    )dnl
])


dnl Check for functioning `signed' keyword
dnl
AC_DEFUN([UC_SIGNED], [dnl
    AC_MSG_CHECKING(for C signed)
    AC_TRY_COMPILE(
        , 
        [signed char x;],
        [AC_MSG_RESULT(yes)],
        [
            AC_MSG_RESULT(no)
            UC_REPLACE(UD_NO_SIGNED,1)
        ]
    )dnl
])


dnl Check for function prototypes
dnl
AC_DEFUN([UC_PROTOTYPES], [dnl
    AC_MSG_CHECKING(for C function prototypes)
    AC_TRY_COMPILE(
        , 
        [extern int foo(int bar);],
        [AC_MSG_RESULT(yes)],
        [
            AC_MSG_RESULT(no)
            UC_REPLACE(UD_NO_PROTOTYPES,1)
        ]
    )dnl
])


dnl Check for the C position-independent-code compile option for building
dnl sharable libraries.
dnl
AC_DEFUN([UC_CPICOPT], [dnl
    AC_REQUIRE([UC_OS]) dnl
    AC_MSG_CHECKING(for C position-independent-code compile-option)
    case "${CPICOPT-unset}" in
        unset)
            case "$GCC" in
                1)  CPICOPT=-fpic;;
                *)  case "$OPSYS" in
                        hpux*)  CPICOPT=+z;;
                        sunos*)
                            case "$OS_MAJOR" in
                                4*)     CPICOPT=-pic;;
                                5*)     CPICOPT=-Kpic;;
                            esac;;
                        *)      CPICOPT=;;
                    esac
                    ;;
            esac
            orig_cflags=$CFLAGS
            CFLAGS="$CFLAGS $CPICOPT"
            AC_TRY_COMPILE(, , , CPICOPT=)
            CFLAGS=$orig_cflags
            ;;
    esac
    AC_MSG_RESULT($CPICOPT)
    AC_SUBST(CPICOPT)
])


dnl Check for the FORTRAN position-independent-code compile option for 
dnl building sharable libraries.
dnl
AC_DEFUN([UC_FPICOPT], [dnl
    AC_REQUIRE([UC_OS]) dnl
    AC_MSG_CHECKING(for FORTRAN position-independent-code compile-option)
    case "${FPICOPT-unset}" in
        unset)
            case "$OPSYS" in
                hpux*)  FPICOPT=+z;;
                sunos*)
                    case "$OS_MAJOR" in
                        4*)     FPICOPT=-pic;;
                        5*)     FPICOPT=-Kpic;;
                    esac;;
                *)      FPICOPT=;;
            esac
            ;;
    esac
    AC_MSG_RESULT($FPICOPT)
    AC_SUBST(FPICOPT)
])


dnl Check for the position-independent-code compile options for 
dnl building sharable libraries.
dnl
AC_DEFUN([UC_PICOPT], [dnl
    AC_REQUIRE([UC_CPICOPT])dnl
    AC_REQUIRE([UC_FPICOPT])dnl
])


dnl Convert argument to uppercase.
dnl
AC_DEFUN([UC_UPPERCASE],[translit($1,abcdefghijklmnopqrstuvwxyz,ABCDEFGHIJKLMNOPQRSTUVWXYZ)])


dnl Return the C macro name version of the argument.
dnl
AC_DEFUN([UC_C_MACRONAME], [UC_UPPERCASE([translit($1,/.<>,__)])])


dnl Check for cpp(1).  This macro replaces the ac_prog_cpp macro because:
dnl     1. That macro, for some reason, sets the value of the shell 
dnl        variable `CPP' to `${CC-cc} -E' rather than to the cpp(1)
dnl        program and such a value has caused trouble in shell command
dnl        lines;
dnl     2. The documentation explicitly states that the ac_prog_cpp macro 
dnl        should be called after the ac_prog_cc macro, so there's no reason 
dnl        for the above value that I can see; and
dnl     3. We need to discover when ${CPP} doesn't work (e.g. when it's 
dnl        defined as `acc -E' under older versions of SunOS).
dnl
AC_DEFUN([UC_PROG_CPP],
[dnl
    AC_REQUIRE([UC_PROG_CC])dnl
    AC_PROG_CPP dnl
    AC_MSG_CHECKING(the C preprocessor)
    AC_TRY_CPP([#include <stdlib.h>],
               [AC_MSG_RESULT(works)],
               [
                   AC_MSG_WARN([[C preprocessor, \`$CPP', doesn't work]])
                   UC_NEED_VALUE(CPP, [C preprocessor], /lib/cpp)
               ])
    AC_SUBST(CPP)
])


dnl Obtain the pathname of a system-supplied header file.  The value of the
dnl associated C macro is `/dev/null' if the header-file could not be
dnl found.
dnl NB: We don't ac_require(uc_prog_cpp) here because it's too deep.
dnl
AC_DEFUN([UC_SYSTEM_H_PATH], [dnl
    echo "#include <$1.h>" > conftestpath.c
    #
    # We add additional `/'s to the header file name to preclude compiler 
    # warnings about the non-portability of `#include "/usr/include/..."'.
    #
    case `uname -s`${GCC-} in
    AIX)
        #
        # AIX's C compiler doesn't emit a line that gives the pathname of
        # the included file.
        #
        # AIX's C compiler puts dependency information in a `.u' file.
        #
        # AIX 4.1's cc(1) makes the following complaint:
        #     ld: 0711-317 ERROR: Undefined symbol: .main
        #
        # AIX 4.1's ksh(1) has problems with `2>&5' so we redirect to 
        # /dev/null.
        #
        $CC_MAKEDEPEND $CPPFLAGS conftestpath.c 2>/dev/null
        path=`sed -n '/[[^\/]]*\(\/[[^  ]]*$1\.h\).*/s//\1/p' \
                conftestpath.u | head -1`
        rm conftestpath.u
        ;;
    HP-UX)
        #
        # HP-UX's C compiler doesn't have a dependency-generation option,
        # so we use another method.
        #
        path=`$CC -E $CPPFLAGS conftestpath.c 2>&5 | 
              sed -n '/[[^\/]]*\(\/[[^  ]]*$1\.h\).*/s//\1/p' |
              head -1`
        ;;
    *)
        path=`$CC_MAKEDEPEND $CPPFLAGS conftestpath.c 2>&5 |
              sed -n '/[[^\/]]*\(\/[[^  ]]*$1\.h\).*/s//\1/p' |
              head -1`
        ;;
    esac
    case "$path" in
    '')
        path=/dev/null
        ;;
    *)
        path=//$path
        ;;
    esac
    AC_DEFINE_UNQUOTED([UD_SYSTEM_[]UC_C_MACRONAME(ifelse($2,,$1,$2))_H],
                       "$path")dnl
])


dnl Obtain the pathname of a system-supplied header file and create a
dnl header file in the port/misc subdirectory that references it.
dnl
AC_DEFUN([UC_SYSTEM_HEADER],
[dnl
    ifelse(
        $2,
        ,
        [UC_SYSTEM_H_PATH($1)],
        [UC_SYSTEM_H_PATH($1, $2)]
    ) dnl
    ifelse(
        $1,
        float,
        ,
        [
            UC_CONFIG_HEADER(port/misc/ifelse($2,,$1,$2).h)
            UC_ENSURE(PORT_SUBDIRS, misc)dnl
        ] dnl
    )dnl
])


dnl Define macros for variadic function support
dnl
AC_DEFUN([UC_VARIADIC_FUNCTIONS],[dnl
    UC_ENSURE(PORT_MANIFEST, stdarg.h.in)dnl
    AC_MSG_CHECKING(for standard C variadic functions)
    AC_TRY_COMPILE([#include <stdarg.h>],
        [
            }
            int foo(int bar, ...) {
            va_list     alist;
            va_start(alist, bar);
            bar = (int)va_arg(alist, int);
            va_end(alist);
        ],
        [
            AC_MSG_RESULT(yes)
            UC_SYSTEM_HEADER(stdarg)
        ],
        [
            AC_MSG_RESULT(no)
            AC_DEFINE(UD_NO_STDARG)dnl
            UC_SYSTEM_HEADER(varargs, stdarg)
        ]
    )dnl
])


dnl Define macro for string generation
dnl
AC_DEFUN([UC_MAKESTRING], [dnl
    AC_MSG_CHECKING(for standard C string generation)
    AC_TRY_COMPILE([# define MAKESTRING(x)      #x],
        char *cp = MAKESTRING(foo);,
        AC_MSG_RESULT(yes),
        [
            AC_MSG_RESULT(no)
            UC_REPLACE(UD_NO_HASH,1)
        ]
    )dnl
])


dnl Define macro for token pasting.
dnl
AC_DEFUN([UC_GLUE], [dnl
    AC_MSG_CHECKING(for standard C token pasting)
    AC_TRY_COMPILE([#define GLUE(a,b) a ## b],
        char *GLUE(c,p) = "foo";,
        AC_MSG_RESULT(yes),
        [
            AC_MSG_RESULT(no)
            UC_REPLACE(UD_NO_HASHHASH,1)
        ]
    )
])


dnl Define pointer-to-void macro.
dnl
AC_DEFUN([UC_VOIDP], [dnl
    AC_MSG_CHECKING(C void pointer)
    AC_TRY_COMPILE(, 
        extern void *foo();,
        AC_MSG_RESULT(yes),
        [
            AC_MSG_RESULT(no)
            UC_REPLACE(UD_NO_VOIDSTAR,1)
        ]
    )
])


dnl pure CFORTRAN support:
dnl
AC_DEFUN([UC_CFORTRAN_PURE], [dnl
UC_MAKEFILE(port/cfortran/Makefile) dnl
AC_REQUIRE([UC_GLUE])dnl
AC_REQUIRE([UC_OS])dnl
AC_MSG_CHECKING(style for cfortran.h)
UC_ENSURE(PORT_SUBDIRS, cfortran)dnl
UC_ENSURE(PORT_ALL, cfortran/all)dnl
UC_ENSURE(PORT_INSTALL, cfortran/install)dnl
UC_ENSURE(PORT_CLEAN, cfortran/clean)dnl
UC_ENSURE(PORT_DISTCLEAN, cfortran/distclean)dnl
case "${UD_NO_HASHHASH-}" in
  1)    CFORTRAN_TYPE=Reiser
        ;;
  *)    CFORTRAN_TYPE=Standard
        ;;
esac
AC_SUBST(CFORTRAN_TYPE)
AC_MSG_RESULT($CFORTRAN_TYPE)
])


dnl CFORTRAN support:
dnl
AC_DEFUN([UC_CFORTRAN], [dnl
UC_CFORTRAN_PURE
case "$OPSYS" in
    aix*)
        # Ensure adherence to the trailing-underscore convention.
        UC_ENSURE(CPPFLAGS, -Dextname)dnl
        UC_ENSURE(FFLAGS, -qextname)dnl
        ;;
    hp*)
        # Ensure adherence to the trailing-underscore convention.
        UC_ENSURE(CPPFLAGS, -Dextname)dnl
        UC_ENSURE(FFLAGS, +ppu)dnl
        ;;
esac
])


dnl Check for standard, udposix(3) stuff.
dnl
AC_DEFUN([UC_UDPOSIX],
[dnl
    AC_REQUIRE([UC_PROG_CC])
    AC_REQUIRE([UC_PROG_CC_MAKEDEPEND])
    AC_REQUIRE([UC_PROG_CPP])
    AC_REQUIRE([UC_CONST])
    AC_REQUIRE([UC_VOLATILE])
    AC_REQUIRE([UC_SIGNED])
    AC_REQUIRE([UC_PROTOTYPES])
    AC_REQUIRE([UC_VARIADIC_FUNCTIONS])
    AC_REQUIRE([UC_MAKESTRING])
    AC_REQUIRE([UC_GLUE])
    AC_REQUIRE([UC_VOIDP])
    UC_CONFIG_HEADER(port/misc/udposix.h)
    UC_ENSURE(PORT_MANIFEST, udposix.h.in)
])


dnl Check for a function declaration in a header file.
dnl
AC_DEFUN([UC_FUNC_DECL], [dnl
    AC_CHECK_HEADER(translit($1, <>),
        [dnl Header-file found:
            AC_MSG_CHECKING(C header file $1 for function $2())
            AC_TRY_COMPILE([#include $1
                    extern struct {int a; int b;} *$2();
                ],
                ,
                [dnl
                    AC_MSG_RESULT(undeclared)
                    UC_REPLACE([UD_NO_[]UC_UPPERCASE($2)_DECL],1)
                ],
                [AC_MSG_RESULT(declared)]
            )dnl
        ], dnl Header-file not found:
        [
            UC_REPLACE([UD_NO_[]UC_UPPERCASE($2)_DECL],1)
        ]) dnl
])


dnl Check for a function.
dnl
AC_DEFUN([UC_FUNC], [dnl
AC_REPLACE_FUNCS($2)dnl
UC_FUNC_DECL($1, $2)dnl
])


dnl Check for a type definition.
dnl The following line, which was in the previous version of this macro:
dnl         typedef void $2;
dnl didn't cause the ip26-irix64_6 cc compiler to error-exit.
dnl So we next tried something else (typedef/variable namespace collision):
dnl         int $2;
dnl Unfortunately, this didn't cause the VAX ULTRIX 4.4 /bin/cc to
dnl error-exit.  So now we try it only if grep can't find the typedef.
dnl
AC_DEFUN([UC_TYPEDEF], [dnl
    AC_REQUIRE([UC_PROG_CPP]) dnl
    AC_CHECK_HEADER(translit($1, <>),
        found=yes,
        found=no)
    case $found in
        yes)
            AC_MSG_CHECKING(C header file $1 for typedef $2)
            echo '#include $1' >conftest.c
            if ($CPP conftest.c | grep $2 >/dev/null) 2>&5; then
                AC_MSG_RESULT(declared)
                rm conftest.c
            else
                AC_TRY_COMPILE([#include $1
                        int $2;
                    ],
                    ,
                    [
                        AC_MSG_RESULT(undeclared)
                        AC_DEFINE([UD_NO_[]UC_C_MACRONAME($1[]_$2)])
                    ],
                    [AC_MSG_RESULT(declared)]
                )
            fi
            ;;
        no)
            AC_DEFINE([UD_NO_[]UC_C_MACRONAME($1[]_$2)])
            ;;
    esac
])


dnl Check for a structure definition.
dnl
AC_DEFUN([UC_STRUCT], [dnl
    AC_MSG_CHECKING(C header file $1 for structure $2)
    AC_TRY_COMPILE([#include $1
            struct $2 {char *foo;};
        ],
        ,
        [
            AC_MSG_RESULT(undeclared)
            AC_DEFINE([UD_NO_[]UC_UPPERCASE($2)[]_STRUCT])
        ],
        [AC_MSG_RESULT(declared)]
    )dnl
])


dnl Ensure a macro definition.
dnl
AC_DEFUN([UC_MACRO], [dnl
    AC_MSG_CHECKING(C header file $1 for macro $2)
    AC_TRY_COMPILE(
        [
#include $1
#ifdef $2
      error
#endif
        ],
        ,
        [AC_MSG_RESULT(undefined)],
        [AC_MSG_RESULT(defined)]
    ) dnl
])


dnl Check for a library.  It would have been nice to see if a compile-link-
dnl execute sequence would have worked (via AC_TRY_RUN) but, with dynamic
dnl libraries under SunOS, the link and execution fail due to unresolved 
dnl references.  Ergo, we just check for the regular `.a' file.
dnl
AC_DEFUN([UC_TEST_LIB], [dnl
if test -z "$$1"; then
  for dir in $2; do
    for name in $3; do
      if test -r $dir/lib$name.a; then
        $1="-L$dir -l$name"
        break 2
      fi
    done
  done
  if test -z "$$1"; then
      UC_NEED_VALUE($1, [$4 library], $5)dnl
  fi
fi
AC_SUBST($1)dnl
])


dnl Set up for using udport/ulog.*
dnl
AC_DEFUN([UC_PORT_ULOG], [dnl
UC_ENSURE(PORT_MANIFEST, ulog.c ulog.h inetutil.c inetutil.h)dnl
])


dnl Check for Unidata threads.
dnl
AC_DEFUN([UC_UTHREADS], [dnl
AC_REQUIRE([UC_UDPOSIX_SIGNAL])dnl
AC_REQUIRE([UC_UDPOSIX_TIME])dnl
UC_ENSURE(PORT_SUBDIRS, misc)dnl
UC_ENSURE(PORT_SUBDIRS, pthreads)dnl
UC_ENSURE(PORT_MANIFEST, uthread.c uthread.h.in nanosleep.c) dnl
UC_MAKEFILE(port/pthreads/Makefile) dnl
UC_MAKEFILE(port/pthreads/include/Makefile) dnl
UC_MAKEFILE(port/pthreads/include/pthread/Makefile) dnl
UC_MAKEFILE(port/pthreads/lib/Makefile) dnl
UC_MAKEFILE(port/pthreads/src/Makefile) dnl
UC_MAKEFILE(port/pthreads/src/gmalloc/Makefile) dnl
UC_CONFIG_HEADER(port/misc/uthread.h) dnl
CPP_THREADS=
case `uname -s` in
    SunOS)
        case `uname -r` in
            4*)
                UC_ENSURE(PORT_ALL, misc/all)
                UC_ENSURE(PORT_INSTALL, misc/install)
                UC_ENSURE(PORT_CLEAN, misc/clean)
                UC_ENSURE(PORT_DISTCLEAN, misc/distclean)
                UC_ENSURE(PORT_ALL, pthreads/all)
                UC_ENSURE(PORT_INSTALL, pthreads/install)
                UC_ENSURE(PORT_CLEAN, pthreads/clean)
                UC_ENSURE(PORT_DISTCLEAN, pthreads/distclean)
                UC_CONFIG_HEADER(port/pthreads/include/pthread/errno.h) dnl
                UC_CONFIG_HEADER(port/pthreads/include/pthread/limits.h) dnl
                UC_CONFIG_HEADER(port/pthreads/include/pthread/signal.h) dnl
                UC_CONFIG_HEADER(port/pthreads/include/pthread/unistd.h) dnl
                UC_SYSTEM_H_PATH(errno)dnl
                UC_SYSTEM_H_PATH(limits)dnl
                UC_SYSTEM_H_PATH(unistd)dnl
                UC_DEFAULT(LD_UTHREADS, dnl
                           [-L[]UC_ABSPATH(port/pthreads/lib) -lpthreads])dnl
                AC_DEFINE(UD_THREAD_FSU_6) dnl
                CPP_THREADS=-I[]UC_ABSPATH(port/pthreads/include)
                UD_OS_SUNOS_4=1
                ;;
            5*)
                AC_DEFINE(UD_THREAD_SUNOS_5) dnl
                UC_STRUCT(<time.h>, timespec)dnl

                AC_MSG_CHECKING(C header files <time.h> and <sys/time.h> for structure timespec)
                AC_TRY_COMPILE([
#                                   include <time.h>
#                                   include <sys/time.h>
                                    struct timespec {char *foo;};
                               ],
                               ,
                               [dnl
                                   AC_MSG_RESULT(undeclared)
                                   AC_DEFINE(UD_NO_TIMESPEC_STRUCT)
                               ],
                               [AC_MSG_RESULT(declared)])dnl
                UC_FUNC(<time.h>, nanosleep)dnl
                UC_ENSURE(UC_LIBOBJS, uthread.o)dnl
                UC_ENSURE(PORT_SUBDIRS, misc)dnl
                UC_DEFAULT(LD_UTHREADS, -lthread)dnl
                UD_OS_SUNOS_5=1
                ;;
        esac
        ;;
    ULTRIX)
        AC_DEFINE(UD_THREAD_DEC_4) dnl
        UC_FUNC(<time.h>, nanosleep)dnl
        UC_ENSURE(UC_LIBOBJS, uthread.o)dnl
        UC_ENSURE(PORT_SUBDIRS, misc)dnl
        UC_DEFAULT(LD_UTHREADS, -lcma -li)dnl
        ;;
    *)
        AC_MSG_WARN(Unidata threads not possible on this platform) 
        ;;
esac
AC_SUBST(CPP_THREADS)dnl
AC_SUBST(LD_UTHREADS)dnl
dnl We don't add LD_UTHREADS to LIBS
dnl to allow the user to specify a dummy value without interfering with
dnl compilation and link tests.
])


dnl Ensure a POSIX <limits.h>.
dnl
AC_DEFUN([UC_UDPOSIX_LIMITS], [dnl
    AC_REQUIRE([UC_UDPOSIX])dnl
    UC_ENSURE(PORT_MANIFEST, limits.h.in)dnl
    UC_ENSURE(PORT_SUBDIRS, misc)dnl
    UC_CONFIG_HEADER(port/misc/limits.h) dnl
    UC_SYSTEM_HEADER(limits) dnl
])


dnl Ensure a POSIX <float.h>.  Test for both existance of <float.h> and
dnl that it's complete.
dnl
AC_DEFUN([UC_UDPOSIX_FLOAT], [dnl
    AC_REQUIRE([UC_UDPOSIX])dnl
    UC_ENSURE(PORT_MANIFEST, config.c)dnl
    AC_HEADER_CHECK(
        float.h,
        [dnl
            AC_TEST_CPP(
                [
#include <float.h>
#ifdef DBL_DIG
#include "DBL_DIG is defined"   /* hard error in case warnings suppressed*/
#endif
                ],
                [
                    UC_ENSURE(PORT_HEADERS, float.h) dnl
                    UC_ENSURE(PORT_SUBDIRS, misc)dnl
                    AC_MSG_RESULT(supplying own)
                ]
            )
        ],
        [
            UC_ENSURE(PORT_HEADERS, float.h) dnl
            UC_ENSURE(PORT_SUBDIRS, misc)dnl
            AC_MSG_RESULT(supplying own)
        ]
    )
])


dnl Ensure a POSIX <search.h> (with tsearch() support)
dnl
AC_DEFUN([UC_UDPOSIX_SEARCH], dnl
[
    AC_REQUIRE([UC_UDPOSIX]) dnl
    UC_ENSURE(PORT_MANIFEST, search.h.in tsearch.c tfind.c tdelete.c twalk.c \
        search-node.h) dnl
    UC_CONFIG_HEADER(port/misc/search.h) dnl
    UC_ENSURE(PORT_SUBDIRS, misc)dnl
    UC_SYSTEM_HEADER(search) dnl
    UC_TYPEDEF(<search.h>, ENTRY) dnl
    UC_TYPEDEF(<search.h>, ACTION) dnl
    UC_TYPEDEF(<search.h>, VISIT) dnl
    UC_FUNC(<search.h>, tsearch)dnl
    UC_FUNC(<search.h>, tfind)dnl
    UC_FUNC(<search.h>, tdelete)dnl
    UC_FUNC(<search.h>, twalk)dnl
])


dnl Ensure a POSIX <stdarg.h>.
dnl
AC_DEFUN([UC_UDPOSIX_STDARG], [dnl
AC_REQUIRE([UC_VARIADIC_FUNCTIONS])dnl
])


dnl Ensure a POSIX <stddef.h>.
dnl
AC_DEFUN([UC_UDPOSIX_STDDEF], [dnl
AC_REQUIRE([UC_UDPOSIX]) dnl
UC_ENSURE(PORT_MANIFEST, stddef.h.in) dnl
UC_CONFIG_HEADER(port/misc/stddef.h) dnl
UC_ENSURE(PORT_SUBDIRS, misc)dnl
UC_SYSTEM_HEADER(stddef) dnl
UC_TYPEDEF(<stddef.h>, size_t) dnl
UC_TYPEDEF(<stddef.h>, ptrdiff_t) dnl
])


dnl Ensure a POSIX <stdio.h>.
dnl
dnl This shouldn't be necessary but is because at least one environment
dnl (SunOS 4.1.3 with the bundled compiler or gcc) doesn't declare
dnl fgetc().
dnl
AC_DEFUN([UC_UDPOSIX_STDIO], [dnl
AC_REQUIRE([UC_UDPOSIX])dnl
UC_ENSURE(PORT_MANIFEST, stdio.h.in)dnl
UC_ENSURE(PORT_SUBDIRS, misc)dnl
UC_SYSTEM_HEADER(stdio)dnl
UC_FUNC(<stdio.h>, fgetc)dnl
])


dnl Ensure a POSIX <stdlib.h>.
dnl
AC_DEFUN([UC_UDPOSIX_STDLIB], [dnl
AC_REQUIRE([UC_UDPOSIX])dnl
UC_ENSURE(PORT_MANIFEST, stdlib.h.in atexit.c)dnl
UC_CONFIG_HEADER(port/misc/stdlib.h) dnl
UC_ENSURE(PORT_SUBDIRS, misc)dnl
UC_SYSTEM_HEADER(stdlib)dnl
dnl UC_TYPEDEF(<stdlib.h>, div_t, struct div { int quot; int rem; })dnl
dnl UC_TYPEDEF(<stdlib.h>, ldiv_t, struct ldiv { long quot; long rem; })dnl
UC_TYPEDEF(<stdlib.h>, size_t)dnl
UC_FUNC(<stdlib.h>, atexit, int atexit, (void (*fcn)(void)))dnl
UC_FUNC(<stdlib.h>, getenv)dnl
])


dnl Ensure a POSIX <string.h>.
dnl
AC_DEFUN([UC_UDPOSIX_STRING], [dnl
AC_REQUIRE([UC_UDPOSIX])dnl
UC_ENSURE(PORT_MANIFEST, strerror.c strstr.c string.h.in memmove.c)dnl
UC_ENSURE(PORT_SUBDIRS, misc)dnl
UC_CONFIG_HEADER(port/misc/string.h) dnl
UC_SYSTEM_HEADER(string)dnl
UC_TYPEDEF(<string.h>, size_t)dnl
UC_FUNC(<string.h>, strerror, char *strerror, (int errno))dnl
UC_FUNC_DECL(<string.h>, strchr)dnl
UC_FUNC_DECL(<string.h>, strcpy)dnl
UC_FUNC_DECL(<string.h>, strrchr)dnl
UC_FUNC_DECL(<string.h>, strncpy)dnl
UC_FUNC_DECL(<string.h>, strtok)dnl
UC_FUNC(<string.h>, strstr, char *strstr, 
    (const char *cs\, const char *ct))dnl
UC_FUNC(<string.h>, memmove, void *memmove,
    (void *s1, const void *s2, size_t n))dnl
dnl AC_HAVE_FUNCS(bcopy [[index]] rindex)dnl
])


dnl Ensure a POSIX <regex.h>.
dnl
AC_DEFUN([UC_UDPOSIX_REGEX], [dnl
AC_REQUIRE([UC_UDPOSIX])dnl
AC_REQUIRE([UC_OS])dnl
PORT_CPP_LDM=
])


dnl Check for the LDM's [regexp](3) library library.
dnl
AC_DEFUN([UC_LIB_LDMREGEXP], [dnl
    AC_MSG_CHECKING([for LDM [regexp](3) library]) 
    UC_TEST_LIB(
        LD_LDMREGEXP,
        [dnl
            /usr/local/ldm/lib dnl
            /upc/ldm/lib dnl
            UC_ABSPATH($prefix/lib) dnl
            UC_ABSPATH($prefix/../lib) dnl
            UC_ABSPATH($prefix/../ldm/lib) dnl
        ],
        [[regexp]],
        LDM [[regexp]](3),
        -L/usr/local/unidata/lib -lregexp
    ) dnl
    AC_MSG_RESULT(${LD_LDMREGEXP-}) 
])


dnl Ensure an alloca().
dnl
AC_DEFUN([UC_ALLOCA], [dnl
UC_ENSURE(PORT_MANIFEST, alloca.c)dnl
AC_FUNC_ALLOCA
case "${ALLOCA-}" in
alloca.o)
    UC_ENSURE(UC_LIBOBJS, alloca.o)dnl
    UC_ENSURE(PORT_SUBDIRS, misc)dnl
    ;;
esac
])


dnl Ensure a POSIX <time.h>.
dnl
AC_DEFUN([UC_UDPOSIX_TIME], [dnl
AC_REQUIRE([UC_UDPOSIX])dnl
UC_ENSURE(PORT_MANIFEST, difftime.c strftime.c time.h.in)dnl
UC_CONFIG_HEADER(port/misc/time.h) dnl
UC_SYSTEM_HEADER(time) dnl
UC_TYPEDEF(<time.h>, time_t)dnl
UC_TYPEDEF(<time.h>, size_t)dnl
UC_FUNC(<time.h>, difftime)dnl
UC_FUNC(<time.h>, strftime)dnl
])


dnl Ensure a POSIX <signal.h>.
dnl
AC_DEFUN([UC_UDPOSIX_SIGNAL], [dnl
    AC_REQUIRE([UC_UDPOSIX])dnl
    UC_ENSURE(PORT_MANIFEST,
        signal.h.in sigaction.c sigaddset.c dnl
        sigdelset.c sigemptyset.c sigfillset.c sigismember.c sigpending.c dnl
        sigprocmask.c sigsuspend.c dnl
    )dnl
    UC_CONFIG_HEADER(port/misc/signal.h)dnl
    UC_SYSTEM_HEADER(signal)dnl
    UC_TYPEDEF(<signal.h>, sigset_t)dnl
    UC_TYPEDEF(<signal.h>, sig_atomic_t)dnl
    UC_STRUCT(<signal.h>, sigaction)dnl
    UC_FUNC(<signal.h>,
        sigaction,
        int sigaction,
        (int sig\, const struct sigaction *act\, struct sigaction *oact) dnl
    )dnl
    case "$UD_NO_SIGACTION_DECL" in
        1)
            AC_HAVE_FUNCS(sigvec sigblock sigpause sigsetmask sigstack bsdsigp) dnl
            ;;
    esac
])


dnl Ensure a POSIX <unistd.h>.
dnl
AC_DEFUN([UC_UDPOSIX_UNISTD], [dnl
AC_REQUIRE([UC_UDPOSIX])dnl
UC_ENSURE(PORT_MANIFEST, unistd.h.in)dnl
UC_CONFIG_HEADER(/port/misc/unistd.h)dnl
UC_SYSTEM_HEADER(unistd)dnl
UC_TYPEDEF(<unistd.h>, size_t)dnl
UC_FUNC(<unistd.h>, getlogin)dnl
])


dnl Ensure a <select.h>.
dnl
AC_DEFUN([UC_SELECT], [dnl
UC_ENSURE(PORT_MANIFEST, select.h)dnl
UC_SYSTEM_HEADER(sys/select)dnl
])


dnl Check for the tar utility.
dnl
AC_DEFUN([UC_PROG_TAR], [dnl
AC_REQUIRE([UC_OS]) dnl
AC_MSG_CHECKING(for tar flags)
case "$OPSYS" in
    aix*)
        TARFLAGS=-chf
        ;;
    bsdos|irix*)
        TARFLAGS=-cLof
        ;;
    *)
        TARFLAGS=-chof
        ;;
esac
AC_MSG_RESULT($TARFLAGS)
AC_SUBST(TARFLAGS) dnl
])


dnl Check for means to generate makefile dependencies from C compiler.
dnl
AC_DEFUN([UC_PROG_CC_MAKEDEPEND], [dnl
AC_BEFORE([$0], [UC_SYSTEM_H_PATH])dnl
AC_REQUIRE([UC_OS]) dnl
AC_MSG_CHECKING(for dependency generation mechanism) 
if test -z "$CC_MAKEDEPEND"; then
    case ${GCC-}${OPSYS} in
        convexos)
            CC_MAKEDEPEND="$CC -k"
            ;;
        hpux)
            CC_MAKEDEPEND="$CC -c -Wp,-M"
            ;;
        sunos5)
            CC_MAKEDEPEND="$CC -xM"
            ;;
        unicos)
            #
            # UNICOS's C compiler has an unusual way of invoking the
            # dependency-generation option.  Also, the c89(1) compiler
            # doesn't even have way of doing this.
            # 
            #
            CC_MAKEDEPEND="cc -c -Wp,-M"
            ;;
        ultrix)
            CC_MAKEDEPEND="$CC -Em"
            ;;
        *)
            CC_MAKEDEPEND="$CC -M"
            ;;
    esac
fi
AC_SUBST(CC_MAKEDEPEND) dnl
AC_MSG_RESULT($CC_MAKEDEPEND) 
])


dnl Check for which(1)
dnl
AC_DEFUN([UC_PROG_WHICH], [dnl
UC_ENSURE(PORT_MANIFEST, which)dnl
AC_CHECK_PROG(WHICH, which, which, [UC_ABSPATH(port)]/which)dnl
])


dnl Check for FORTRAN compiler.
dnl
AC_DEFUN([UC_PROG_FC], [dnl
if test -z "${FC+set}"; then
  AC_REQUIRE([UC_OS])dnl
  case "$OPSYS$OS_MAJOR" in
    aix*)
        AC_PROGRAMS_CHECK(FC, $fc f77 cf77)dnl
        ;;
    hpux*)
        AC_PROGRAMS_CHECK(FC, $fc fort77 fortc f77)dnl
        UC_ENSURE(FFLAGS, +U77)dnl
        UC_ENSURE(LD_FORTRAN, -lU77)dnl
        ;;
    dgux*)
        AC_PROGRAMS_CHECK(FC, $fc ghf77 f77)dnl
        UC_DEFAULT(LD_FORTRAN)dnl
        ;;
    sunos*)
        AC_PROGRAMS_CHECK(FC, $fc f77 cf77)dnl
        UC_DEFAULT(LD_FORTRAN)dnl
        ;;
    irix5*)
        AC_PROGRAMS_CHECK(FC, $fc f77 cf77)dnl
        UC_DEFAULT(LD_FORTRAN, -tl -Bstatic)dnl
        ;;
    *)
        AC_PROGRAMS_CHECK(FC, $fc f77 cf77 fort77)dnl
        UC_DEFAULT(LD_FORTRAN)dnl
        ;;
  esac
  if test -z "$FC"; then
    UC_NEED_VALUE(FC, [FORTRAN compiler], /bin/f77)dnl
  fi
else
  AC_MSG_CHECKING(for FORTRAN compiler)
  AC_MSG_RESULT($FC)
fi
AC_SUBST(FC)
])


dnl Check for FORTRAN library.
dnl
AC_DEFUN([UC_LIB_F77], [dnl
AC_REQUIRE([UC_PROG_FC])dnl
AC_REQUIRE([UC_PROG_WHICH])dnl
AC_MSG_CHECKING(for FORTRAN library) 
case `$WHICH "$FC"` in
  *lang*)
        LD_F77='-lF77 -lM77';;
  *)
        LD_F77='-lF77';;
esac
AC_SUBST(LD_F77)dnl
AC_MSG_RESULT($LD_F77) 
])


dnl Check size of FORTRAN REAL.
dnl
AC_DEFUN([UC_CHECK_SIZEOF_REAL], [dnl
    case "$SIZEOF_REAL" in
    '')
        AC_REQUIRE([UC_PROG_FC])dnl
        AC_REQUIRE([UC_PROG_CC])dnl
        AC_MSG_CHECKING(size of FORTRAN REAL)
        case "$FC" in
        '') ;;
        *)
            cat >conftest.c <<EOF
#               include <stdio.h>
                void sub(p1,p2) char *p1, *p2;
                    { (void)printf("%ld\n", (long)(p2-p1)); }
                void sub_(p1,p2) char *p1,*p2; { sub(p1,p2); }
                void SUB(p1,p2) char *p1,*p2; { sub(p1,p2); }
                void SUB_(p1,p2) char *p1,*p2; { sub(p1,p2); }
EOF
            if $CC -c $CFLAGS $CPPFLAGS conftest.c 2>&5; then
                mv conftest.o conftest.c.o
                cat >conftest.f <<EOF
                    REAL VEC(2)
                    CALL SUB(VEC(1), VEC(2))
                    END
EOF
                if $FC -o conftest $FFLAGS $LDFLAGS conftest.f conftest.c.o \
                      2>&5; then
                    SIZEOF_REAL=`./conftest`
                    break;
                fi
            fi
            AC_SUBST(SIZEOF_REAL)
            rm -rf conftest conftest.o conftest.c.o conftest.f conftest.c
            ;;
        esac
        ;;
    esac
    AC_MSG_RESULT($SIZEOF_REAL)
])


dnl Check for yacc(1) utility.
dnl
AC_DEFUN([UC_PROG_YACC], [dnl
AC_PROGRAMS_CHECK(YACC, yacc bison)
case "${YACC-}" in
    '')
        UC_NEED_VALUE(YACC, [grammar compiler], /usr/bin/yacc)dnl
        ;;
    *bison*)
        UC_ENSURE(YACC, -y)
        ;;
esac
AC_SUBST(YACC)
])


dnl Check for yacc(1) library.
dnl
AC_DEFUN([UC_LIB_YACC],
[dnl
    case `uname` in
        Linux)
            UC_DEFAULT(LD_YACC, )
            ;;
        *)  UC_CHECK_LIB(LD_YACC, yyerror(""), , y, yacc, -ly)
            ;;
    esac
])


dnl Check for library utility, ar(1).
dnl
AC_DEFUN([UC_PROG_AR], [dnl
# We use ar(1) under UNICOS even though bld(1) is preferred because bld(1)
# doesn't understand the "u" option.
AC_CHECK_PROG(AR, ar, ar, )dnl
if test -z "$AR"; then
  UC_NEED_VALUE(AR, [library utility], /bin/ar)dnl
fi
AC_SUBST(AR)
])


dnl Check for troff(1).
dnl
AC_DEFUN([UC_PROG_TROFF], [dnl
AC_CHECK_PROG(TROFF, troff, ptroff, troff)dnl
if test -z "$TROFF"; then
  UC_NEED_VALUE(TROFF, [troff(1)-like utility], /bin/troff)dnl
fi
AC_SUBST(TROFF)
])


dnl Check for fortc(1)
dnl
AC_DEFUN([UC_PROG_FORTC], [dnl
AC_REQUIRE([UC_OS])dnl
AC_REQUIRE([UC_UDPOSIX_STDDEF])dnl
UC_ENSURE(PORT_SUBDIRS, fortc misc)dnl
UC_ENSURE(PORT_ALL, fortc/all)dnl
UC_ENSURE(PORT_INSTALL, fortc/install)dnl
UC_ENSURE(PORT_CLEAN, fortc/clean)dnl
UC_ENSURE(PORT_DISTCLEAN, fortc/distclean)dnl
UC_ENSURE(PORT_MANIFEST, udalloc.c udalloc.h)dnl
UC_MAKEFILE(port/fortc/Makefile)dnl
dir=`pwd`/port/fortc
FORTC="$dir/fortc"
AC_SUBST(FORTC)dnl
])


dnl Check for neqn(1).
dnl
AC_DEFUN([UC_PROG_NEQN], [dnl
AC_CHECK_PROG(NEQN, neqn, neqn, cat)dnl
test "$NEQN" = cat && 
  AC_MSG_WARN($[]0: Can't find program \`neqn'; setting to \`cat') 
])


dnl Check for tbl(1).
dnl
AC_DEFUN([UC_PROG_TBL], [dnl
AC_CHECK_PROG(TBL, tbl, tbl, cat)dnl
test "$TBL" = cat && 
  AC_MSG_WARN($[]0: Can't find program \`tbl'; setting to \`cat') 
])


dnl Check for makeinfo(1).
dnl
AC_DEFUN([UC_PROG_MAKEINFO], [dnl
AC_CHECK_PROG(MAKEINFO, makeinfo, makeinfo)dnl
])


dnl Determine the machine type.
dnl
AC_DEFUN([UC_MACHINE], [dnl
    AC_REQUIRE([UC_OS])dnl
    AC_MSG_CHECKING(type of machine)
    if test -z "$MACHINE"; then
    MACHINE=`uname -m | tr [[A-Z]] [[a-z]]`
    case $OPSYS in
        aix*)
            MACHINE=rs6000
            ;;
        hp*)
            MACHINE=hp`echo $MACHINE | sed 's%/.*%%'`
            ;;
        sunos*)
            case $MACHINE in
                sun4*)
                    MACHINE=sun4
                    ;;
            esac
            ;;
        irix*)
            case $MACHINE in
                ip20)
                    MACHINE=sgi
                    ;;
            esac
            ;;
    esac
    if test -z "$MACHINE"; then
      UC_NEED_VALUE(MACHINE, [machine hardware type], sun4)dnl
    fi
    fi
    AC_SUBST(MACHINE)dnl
    AC_MSG_RESULT($MACHINE) 
])


dnl Check for ncdump(1)
dnl
AC_DEFUN([UC_PROG_NCDUMP], [dnl
AC_PROGRAMS_CHECK(NCDUMP, ncdump [UC_ABSPATH($exec_prefix)]/bin/ncdump)dnl
if test -z "$NCDUMP"; then
  UC_NEED_VALUE(NCDUMP, [netCDF lister], /usr/local/unidata/bin/ncdump)dnl
fi
AC_SUBST(NCDUMP)
])


dnl Check for ncgen(1)
dnl
AC_DEFUN([UC_PROG_NCGEN], [dnl
AC_PROGRAMS_CHECK(NCGEN, ncgen [UC_ABSPATH($exec_prefix)]/bin/ncgen)dnl
if test -z "$NCGEN"; then
  UC_NEED_VALUE(NCGEN, [netCDF generator], /usr/local/unidata/bin/ncgen)dnl
fi
AC_SUBST(NCGEN)
])


dnl Test a script.
dnl
AC_DEFUN([UC_TEST_SCRIPT], [dnl
cat << UD_EOF > conftest.sh
[$1]
UD_EOF
chmod +x conftest.sh
if ./conftest.sh 2>&5; then
  ifelse([$2], , :, [$2])
ifelse([$3], , , [else
  $3
])dnl
fi
rm -f conftest.sh
])dnl


dnl Filter a file through cpp(1).
dnl
AC_DEFUN([UC_FILTER_CPP], [dnl
    AC_REQUIRE([UC_PROG_CPP])dnl
    AC_MSG_WARN(processing $1 with the C preprocessor to produce $2) 
    ifdef([AC_CONFIG_NAME],
        [
            UC_TEST_SCRIPT(
                [dnl
                    echo "$DEFS" > conftest.c
                    echo "# line 1 $1" >> conftest.c
                    cat $1 >> conftest.c
                    $CPP conftest.c | \
                        awk '/^$/ {if (set) next; set=1} {print} !/^$/ {set=0}'\
                            > $2
                    rm -f conftest.c
                ]
            )
        ],
        [
            UC_TEST_SCRIPT(
                [
                    $CPP "$DEFS" $1 | \
                        awk '/^$/ {if (set) next; set=1} {print} !/^$/ {set=0}'\
                            > $2
                ]
            )
        ]
    )
])


dnl Convert a pathname to an absolute one at autoconf(1) execution time.
dnl
AC_DEFUN([UC_ABSPATH_M4], [dnl
syscmd([case "$1" in 
  /*)
        echo $1; exit;;
   *)
        path=`pwd`/$1
        tail=
        while test "$path"; do
            (cd $path && echo `pwd`$rest) && exit
            base=/`basename "$path"`
            tail=/$base$tail
            path=`echo "$path" | sed "s/\/$base//"`
        done;;
esac > conftest.syscmd 2>&1
])dnl
include(conftest.syscmd)dnl
])


dnl Convert a pathname to an absolute one at ./configure execution time.
dnl NB: The [)] construction is neccessary with GNU m4 1.4.
dnl
AC_DEFUN([UC_ABSPATH], [`dnl
case "$1" in 
  /*[)] echo $1; exit;;
   *[)] path=\`pwd\`/$1
        tail=
        while test "$path"; do
          (cd $path && echo \`pwd\`$rest) && exit
          base=/\`basename "$path"\`
          tail=/$base$tail
          path=\`echo "$path" | sed "s/\/$base//"\`
        done;;
esac
`])


dnl Set a value for the installation prefix.
dnl
AC_DEFUN([UC_PREFIX], [dnl
    AC_BEFORE([$0],[UC_PROG_FORTC])dnl
    AC_BEFORE([$0],[UC_LIB_NETCDF])dnl
    AC_BEFORE([$0],[UC_CPP_NETCDF])dnl
    AC_BEFORE([$0],[UC_LIB_NCOPERS])dnl
    AC_BEFORE([$0],[UC_CPP_NCOPERS])dnl
    AC_BEFORE([$0],[UC_LIB_UDPORT])dnl
    AC_MSG_CHECKING(the installation prefix)
    case ${prefix-} in
        NONE|'') prefix=UC_ABSPATH($1);;
    esac
    AC_MSG_RESULT($prefix)
    AC_MSG_CHECKING(the installation exec-prefix)
    case ${exec_prefix-} in
        NONE|'')
            exec_prefix=$prefix;;
    esac
    AC_MSG_RESULT($exec_prefix)
])


dnl Check for a directory containing a file.
dnl
AC_DEFUN([UC_TEST_DIR], [dnl
  if test -z "$$1"; then
    for dir in $2; do
      if test -r $dir/$3; then
        $1=$dir
        break;
      fi
    done
    if test -z "$$1"; then
      UC_NEED_VALUE($1, $4, $5)dnl
    fi
  fi
AC_SUBST($1)dnl
])


dnl Check for McIDAS header-file directory.
dnl
AC_DEFUN([UC_CPP_MCIDAS], [dnl
AC_MSG_CHECKING(for McIDAS header-files) 
UC_TEST_DIR(CPP_MCIDAS, ~mcidas/include ~mcidas/src /usr/mcidas/include /usr/mcidas/src /usr/local/mcidas/include /usr/local/mcidas/src /usr/mcidas/src /home/mcidas/src /home/mcidasd/src, hex80.inc,
    McIDAS headers, -I/home/mcidas/src)dnl
CPP_MCIDAS=`case ${CPP_MCIDAS} in
                -I*)
                    echo ${CPP_MCIDAS};;
                *)
                    echo -I${CPP_MCIDAS-};;
            esac`
AC_MSG_RESULT($CPP_MCIDAS) 
])


dnl Check for McIDAS library.
dnl
AC_DEFUN([UC_LIB_MCIDAS], [dnl
AC_REQUIRE([UC_USLEEP])dnl
AC_MSG_CHECKING(for McIDAS library) 
UC_TEST_LIB(LD_MCIDAS, ~mcidas/lib /usr/mcidas/lib /usr/local/mcidas/lib /home/mcidas/lib /home/mcidasd/lib, mcidas, McIDAS, dnl
  -L/home/mcidas/lib -lmcidas)dnl
AC_MSG_RESULT($LD_MCIDAS) 
])


dnl Check for X11 header-file directory.
dnl
AC_DEFUN([UC_CPP_X11], [dnl
AC_MSG_CHECKING(for X11 header-files) 
UC_TEST_DIR(CPP_X11, ${OPENWINHOME-/usr/openwin}/[[include]] dnl
    /usr/[[include]] /usr/[[include]] /usr/local/[[include]], X11/Xlib.h,
    X11 header, -I/usr/openwin/[[[include]]])dnl
CPP_X11=`case ${CPP_X11} in
            -I*)
                echo ${CPP_X11};;
            *)
                echo -I${CPP_X11-};;
        esac`
AC_MSG_RESULT($CPP_X11) 
])


dnl Form a library reference for the linker/loader
dnl
dnl On a SunOS 5 system, a `-R<dir>' is added in addition to a `-L<dir>'
dnl in order to make the utility independent of LD_LIBRARY_PATH (is this
dnl a good idea?) and to ensure that it'll run regardless of who
dnl executes it.
dnl
dnl UC_LINK_REF(varname, libdir, libname)
dnl
dnl Example: UC_LINK_REF(UC_LD_MATH, /upc/netcdf/lib, netcdf)
dnl
AC_DEFUN([UC_LINK_REF], [dnl
    AC_REQUIRE([UC_OS]) dnl
    case "${OPSYS}$OS_MAJOR" in
        unicos*)
            case "$2" in
                '') $1="-l $3";;
                *)  $1="-L $2 -l $3";;
            esac
            ;;
        sunos5*)
            case "$2" in
                '') $1="-l$3";;
                *)  $1="-R$2 -L$2 -l$3";;
            esac
            ;;
        *)
            case "$2" in
                '') $1="-l$3";;
                *)  $1="-L$2 -l$3";;
            esac
            ;;
    esac
])


dnl Check for a library that contains a function.
dnl
dnl NB: Always checks default library and library directories first.  This
dnl obviates the need for a `-L...' reference, which can cause problems
dnl (e.g. a `-L/usr/lib -lsocket' reference under SunOS 5.2 can cause the
dnl wrong `-lm' to be loaded).
dnl
dnl This rule was changed (for some reason) to return `-lc' if the 
dnl function was in the default library.  This caused problems on
dnl an DecStation ULTRIX system when f77(1) was used to link a FORTRAN
dnl program: the C and FORTRAN libraries had duplicate definitions for
dnl some functions.  Consequently, we return to the practice of not
dnl deciding on `-lc'.
dnl
dnl UC_CHECK_LIB(varname, func, dir ..., lib ..., libname, example)
dnl
AC_DEFUN([UC_CHECK_LIB],
[dnl
    AC_MSG_CHECKING(for $5 library)
    case "${$1+set}" in
    set)
        AC_MSG_RESULT($$1)
        ;;
    *)
        AC_MSG_RESULT()
        LIBS_save=$LIBS
        found=no
        AC_MSG_CHECKING(for $2 in default library(s))
        AC_TRY_LINK(, $2;, 
                    [
                        AC_MSG_RESULT(yes)
                        $1=
                        found=yes
                    ],
                    [
                        AC_MSG_RESULT(no)
                        os=`uname -s`
                        for dir in '' $3; do
                            for lib in $4; do
                                UC_LINK_REF(LIBS, $dir, $lib)
                                AC_MSG_CHECKING(for $2 in $LIBS)
                                AC_TRY_LINK(, $2;, 
                                            [
                                                AC_MSG_RESULT(yes)
                                                $1=$LIBS
                                                found=yes
                                                break 2
                                            ])
                                AC_MSG_RESULT(no)
                            done
                        done
                    ])
        LIBS=$LIBS_save
        case $found in
            no) UC_NEED_VALUE($1, [$5 library], $6) dnl
                ;;
        esac
        ;;
    esac
    AC_SUBST($1) dnl
])


dnl Check for X11 library.
dnl
AC_DEFUN([UC_LIB_X11],
[dnl
    UC_CHECK_LIB(LD_X11,
                 XNoOp(0),
                 ${OPENWINHOME-/usr/openwin}/lib dnl
                 /usr/lib/X11 /usr/X11/lib /usr/lib/X11R5 dnl
                 /usr/local/lib /usr/local/lib/X11 /usr/local/lib/X11R5,
                 X11,
                 X11,
                 -L/usr/lib/X11 -lX11)dnl
])


dnl Check for usleep().
dnl
AC_DEFUN([UC_USLEEP], [dnl
UC_ENSURE(PORT_MANIFEST, usleep.c)dnl
UC_FUNC(<unistd.h>, usleep)dnl
])


dnl Check for X11 implementation (header file and library).
dnl
AC_DEFUN([UC_X11], [AC_REQUIRE([UC_CPP_X11])AC_REQUIRE([UC_LIB_X11])])


dnl Check for rpc library.
dnl
dnl NB: Don't use clnt_create() for the target because, for some reason,
dnl that yields a false negative under ULTRIX.
dnl
AC_DEFUN([UC_LIB_RPC],
[dnl
    AC_REQUIRE([UC_LIB_SOCKET])
    case "$LD_RPC" in
        '') case `uname -sr` in
                'SunOS 5'*)
                    LD_RPC="-R/usr/ucblib -L/usr/ucblib -lrpcsoc -lnsl $LD_SOCKET"
                    ;;
                *)  UC_CHECK_LIB(LD_RPC, clnttcp_create(),
                                 /usr/lib /usr/local/lib, 
                                 sun nsl rpc, RPC, -L/usr/lib -lnsl)
                    ;;
            esac
            ;;
    esac
])


dnl Check for socket library.
dnl
AC_DEFUN([UC_LIB_SOCKET],
[dnl
    UC_CHECK_LIB(LD_SOCKET, socket(0,0,0), /usr/lib /usr/local/lib, 
                 socket, BSD sockets, 
                 -L/usr/lib -lsocket)dnl
])


dnl Check for Standard C math library.
dnl
AC_DEFUN([UC_LIB_MATH],
[dnl
    UC_CHECK_LIB(LD_MATH, return sin(0.0) != 0, , m, C math, -L/opt/SUNWspro/SC3.0/lib -lm)dnl
])


dnl Check for netCDF header-file directory.
dnl
AC_DEFUN([UC_CPP_NETCDF], [dnl
AC_MSG_CHECKING(for netCDF header-file) 
UC_TEST_DIR(CPP_NETCDF, /upc/netcdf/include dnl
    [UC_ABSPATH($prefix/[[include]])], netcdf.h,
    [netCDF header], [-I/usr/local/unidata/[[include]]])dnl
CPP_NETCDF=`case ${CPP_NETCDF} in -I*) echo ${CPP_NETCDF};; *) echo -I${CPP_NETCDF-};; esac`
AC_MSG_RESULT($CPP_NETCDF) 
])


dnl Check for netCDF library.
dnl
AC_DEFUN([UC_LIB_NETCDF], [dnl
AC_MSG_CHECKING(for netCDF library) 
UC_TEST_LIB(LD_NETCDF, /upc/netcdf/lib dnl
    [UC_ABSPATH($prefix/lib)], netcdf,
  netCDF, -L/usr/local/unidata/lib -lnetcdf)dnl
AC_MSG_RESULT($LD_NETCDF) 
])


dnl Check for netCDF implementation (header file and library).
dnl
AC_DEFUN([UC_NETCDF], [AC_REQUIRE([UC_CPP_NETCDF])AC_REQUIRE([UC_LIB_NETCDF])])


dnl Check for netCDF operators library.
dnl
AC_DEFUN([UC_LIB_NCOPERS], [dnl
AC_MSG_CHECKING(for netCDF operators library) 
UC_TEST_LIB(LD_NCOPERS, [UC_ABSPATH($prefix/lib)], ncopers,
  netCDF-operators, [-L/usr/local/unidata/lib -lncopers])dnl
AC_MSG_RESULT($LD_NCOPERS) 
])


dnl Check for udunits header-file directory.
dnl
AC_DEFUN([UC_CPP_UDUNITS], [dnl
AC_MSG_CHECKING(for udunits header-file) 
UC_TEST_DIR(CPP_UDUNITS, /upc/udunits/include dnl
    [UC_ABSPATH($prefix/[[include]])], udunits.h,
    [udunits header], [-I/usr/local/unidata/[[include]]])dnl
CPP_UDUNITS=`case ${CPP_UDUNITS} in -I*) echo ${CPP_UDUNITS};; *) echo -I${CPP_UDUNITS-};; esac`
AC_MSG_RESULT($CPP_UDUNITS) 
])


dnl Check for udunits library.
dnl
AC_DEFUN([UC_LIB_UDUNITS], [dnl
AC_MSG_CHECKING(for udunits library) 
UC_TEST_LIB(LD_UDUNITS,  /upc/udunits/lib dnl
    [UC_ABSPATH($prefix/lib)], udunits,
    udunits, -L/usr/local/unidata/lib -ludunits)dnl
case "$LD_UDUNITS" in
    '') ;;
    *)  LD_UDUNITS="$LD_UDUNITS -ludport";;
esac
AC_MSG_RESULT($LD_UDUNITS) 
])

dnl Check for udunits implementation (header file and library).
dnl
AC_DEFUN([UC_UDUNITS], [AC_REQUIRE([UC_CPP_UDUNITS])AC_REQUIRE([UC_LIB_UDUNITS])])


dnl Check for LDM header-file directory.
dnl
AC_DEFUN([UC_CPP_LDM], 
[dnl
    AC_BEFORE([$0], [UC_UDPOSIX_REGEX]) dnl
    AC_MSG_CHECKING(for LDM header-file) 
    UC_TEST_DIR(CPP_LDM,
                /usr/local/ldm/[[include]] dnl
                    /upc/ldm/[[include]] dnl
                    [UC_ABSPATH($prefix/[[include]])] dnl
                    [UC_ABSPATH($prefix/../[[include]])] dnl
                    [UC_ABSPATH($prefix/../ldm/[[include]])],
                ldm.h,
                [LDM header],
                [-I/usr/local/unidata/[[include]]]) dnl
    CPP_LDM=`case ${CPP_LDM} in
                 -I*) echo ${CPP_LDM};;
                 *)   echo -I${CPP_LDM-};;
             esac`
    AC_MSG_RESULT($CPP_LDM) 
])


dnl Check for LDM library.
dnl
AC_DEFUN([UC_LIB_LDM], [dnl
AC_MSG_CHECKING(for LDM library) 
UC_TEST_LIB(LD_LDM,
            /usr/local/ldm/lib dnl
            /upc/ldm/lib dnl
    [UC_ABSPATH($prefix/lib)] dnl
                [UC_ABSPATH($prefix/../lib)] dnl
                [UC_ABSPATH($prefix/../ldm/lib)],
            ldm,
            LDM,
            -L/usr/local/unidata/lib -lldm)dnl
AC_MSG_RESULT($LD_LDM) 
])


dnl Check for LDM implementation (header file and library).
dnl
AC_DEFUN([UC_LDM], [AC_REQUIRE([UC_CPP_LDM])AC_REQUIRE([UC_LIB_LDM])])


dnl Check for udres(3) library.
dnl
AC_DEFUN([UC_LIB_UDRES], [dnl
AC_MSG_CHECKING(checking for udres library) 
UC_TEST_LIB(LD_UDRES, [UC_ABSPATH($prefix/lib)], udape,
  udres, -L/usr/local/unidata/lib -ludape)dnl
AC_MSG_RESULT($LD_UDRES) 
])


dnl Set installation programs.  This differs from the standard
dnl autoconf(1) macro by making installed data files group writable.
dnl
AC_DEFUN([UC_PROG_INSTALL], [dnl
AC_PROG_INSTALL[]dnl
INSTALL_DATA="`echo "${INSTALL_DATA}" | sed 's/644/664/'`"
])


dnl Check for missing definitions.
dnl
AC_DEFUN([UC_CHECK_MISSING], [dnl
if test -s confdefs.missing; then
  echo
  echo "$[]0: ERROR: The following variables need values:"
  echo
  awk -F: 'BEGIN {printf "%-13s%-27s%s\n", "VARIABLE", "DESCRIPTION", "EXAMPLE";
                  printf "%-13s%-27s%s\n", "--------", "-------", "-------"}
                 {printf "%-13s%-27s%s\n", $[]1, $[]2, $[]3}' confdefs.missing
  if test -t 0 -a -t 1; then
    cat << \UD_CAT_EOF

For each variable above, this script will now request that you input
a value appropriate for your system (note that this value will not
be interpolated by a shell -- so don't use shell substitutions).
Alternatively, you can interrupt this script, set the above variables
in the environment or in the file CUSTOMIZE, and then re-execute this
script.  (Variables referring to executable programs needn't be set if
the relevant directory is added to PATH.  See file INSTALL for details.)

UD_CAT_EOF
    saveifs="$IFS"; IFS=:
    while read variable description example; do
      echo "Enter a value for the $description (e.g. \"$example\"):" >/dev/tty
      read value </dev/tty
      echo "$variable='$value'"
      echo >/dev/tty
    done <confdefs.missing >conf.assigns
    IFS="$saveifs"
    . ./conf.assigns
    rm -f conf.assigns
  else
    cat << UD_CAT_EOF

The above variables may be set in the environment or in the file
CUSTOMIZE.  Variables referring to executable programs needn't be set
if the relevant directory is added to PATH.  In any case, ./configure
should probably be rerun.  See file INSTALL for details.

UD_CAT_EOF
    rm confdefs.missing
    exit 1
  fi
fi
rm -f confdefs.missing
])


dnl Post process makefiles.
dnl
AC_DEFUN([UC_POSTPROCESS_MAKEFILES], 
[dnl
dnl
dnl Don't put an AC_REQUIRE([UC_PROG_CC]) here because the context is
dnl one of being evaluated by AC_OUTPUT() and it messes up.
dnl
#
# Create a script to accomplish the post processing.
#
cat << UD_EOF_CONFTEST_C > conftest.c
changequote(<<,>>)dnl
#include <stdio.h>
int readsub(inpath)
    char        *inpath;
{
    char        buf[2048], path[1024];
    FILE        *fp     = inpath == NULL
                                ? stdin
                                : fopen(inpath, "r");
    if (fp == NULL) {
        (void) perror(inpath);
        return 0;
    }
    buf[sizeof(buf)-1]  = 0;
    while (fgets(buf, sizeof(buf), fp) != NULL) {
        if (sscanf(buf, "include%*[ \t]%s", path) == 1) {
            if (!readsub(path))
                return 0;
        } else {
            (void) fputs(buf, stdout);
        }
    }
    return 1;
}
int main(int argc, char* argv[])
{
    return readsub((char*)NULL) ? 0 : 1;
}
changequote([,])dnl
UD_EOF_CONFTEST_C
if eval $CC -o conftest ${CFLAGS-} ${CPPFLAGS-} ${LDFLAGS-} conftest.c; then
    conftest=`pwd`/conftest
    case "$1" in
    '')
        ;;
    *)
        for file in $1; do
            echo 1>&2 "expanding \`include's in file \`$file'"
            sd=`pwd`/`echo $file | sed 's,[[^/]]*$,,'`
            base=`basename $file`
            (cd $sd && $conftest < $base > conftest.mk && 
              mv conftest.mk $base) || exit 1
        done
        ;;
    esac
fi
rm conftest conftest.c
])


dnl Get shell-variable override values for local customizations.
dnl
AC_DEFUN([UC_CUSTOMIZE],
[dnl
    AC_BEFORE([$0], [UC_DEFAULT])dnl
    if test -r CUSTOMIZE; then
        . ./CUSTOMIZE
    fi
])


dnl Set the root of the FTP distribution directory.
dnl
AC_DEFUN([UC_FTPDIR], [dnl
AC_BEFORE([UC_BINFTP])
FTPDIR=${FTPDIR-/web/ftp/pub}/$1
AC_SUBST(FTPDIR)dnl
])


dnl Set the binary distribution directory.
dnl
AC_DEFUN([UC_BINFTP], [dnl
    AC_MSG_CHECKING([binary distribution directory])
    case ${FTPBINDIR-unset} in
        unset)
            system=`system 2>/dev/null || echo dummy_system`
            FTPBINDIR=${FTPDIR-/home/ftp}/pub/binary/$system
            ;;
    esac
    AC_SUBST(FTPBINDIR)dnl
    AC_MSG_RESULT($FTPBINDIR)
])


dnl Set package name.
dnl
AC_DEFUN([UC_PACKAGE], [dnl
AC_MSG_CHECKING(for package name) 
PACKAGE=${PACKAGE-`basename \`pwd\``}
AC_SUBST(PACKAGE)dnl
AC_MSG_RESULT($PACKAGE) 
])


dnl Set version identifiers.
dnl
AC_DEFUN([UC_VERSION], [dnl
AC_MSG_CHECKING(for package version) 
if test -z "$VERSION"; then
  if test -r VERSION; then \
    VERSION="`cat VERSION`"
  else
    VERSION=
  fi
fi
AC_SUBST(VERSION)dnl
if test -z "$MAJOR_NO"; then
  if test "$VERSION"; then \
    MAJOR_NO=`echo $VERSION |
      sed -n '/^\([[0-9]][[0-9]]*\)\.[[0-9]][[0-9]]*.*/s//\1/p;q'`
  else
    MAJOR_NO=
  fi
fi
AC_SUBST(MAJOR_NO)dnl
if test -z "$MINOR_NO"; then
  if test "$VERSION"; then \
    MINOR_NO=`echo $VERSION |
      sed -n '/^[[0-9]][[0-9]]*\.\([[0-9]][[0-9]]*\).*/s//\1/p;q'`
  else
    MINOR_NO=
  fi
fi
AC_SUBST(MINOR_NO)dnl
AC_MSG_RESULT($MAJOR_NO.$MINOR_NO) 
])


dnl Handle a missing value.
dnl
AC_DEFUN([UC_NEED_VALUE], [dnl
echo "$1:$2:$3" >> confdefs.missing
])


dnl Ensure that a variable contains a given string and that it's substituted.
dnl NB: The [)] construction is neccessary with GNU m4 1.4.
dnl
AC_DEFUN([UC_ENSURE], [dnl
ifelse($2, , [dnl
  $1=${$1-}
], [dnl
  for arg in $2; do
    case "$$1" in
      *$arg*[)]
        ;;
      *[)]
        $1="${$1-} $arg"
        ;;
    esac
  done
])dnl
AC_SUBST($1)dnl
]dnl
)
