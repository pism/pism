.. default-role:: literal

PISM coding guidelines
======================

.. contents::

File names
----------

C++ source files use the extension `.cc`. Headers use `.hh`. C code uses `.c` and `.h`.

The *basename* of a file containing the source code for a class should match the name of
this class, including capitalization. For example, a class `Foo` is declared in `Foo.hh`
and implemented in `Foo.cc`.

Source and header files
-----------------------

Each header file must have an include guard. Do *not* use "randomized" names of include
guards.

Headers should *not* contain function and class method implementations unless these
implementations *should be inlined*; see below.

Inline functions and methods
----------------------------

Implementations of inline methods should be *outside* of class declarations and in a
*separate header* called `Foo_inline.hh`. This header will then be included from `Foo.hh`.

Including header files
----------------------

Include all system headers before PISM headers. Separate system headers from PISM headers
with an empty line.

Good:

.. code-block:: c++

   #include <cassert>
   #include <cstring>

   #include "IceGrid.hh"

Bad:

.. code-block:: c++

   #include <cstring>
   #include "IceGrid.hh"
   #include <cassert>

Whenever appropriate add comments explaining why a particular header was included.

.. code-block:: c++

   #include <cstring> // strcmp

Naming
------

Variable
^^^^^^^^

- Variable names should use lower case letters and (if necessary) digits separated by
  underscores, for example: `ice_thickness`.
- Do not abbreviate names of variables used in more than one scope **unless** this is
  needed to keep the name under 30 characters. If a function uses a variable so many times
  typing a long name is a hassle, create a reference with a shorter name that is only
  visible in this scope. (This alias will be compiled away.)
- Single-letter variable names should only be used in "small" scopes (short functions,
  etc.)
- If you are going to use a single-letter name, pick a letter that is easy to read (avoid
  `i`, `l`, and `o`).
- Names of class data members should use the `m_` prefix, for example: `m_name`.
- Names of static data members should use the `sm_` prefix.
- Global variables (which should be avoided in general) use the `g` prefix.

Types and classes
^^^^^^^^^^^^^^^^^

Names of types and classes use `CamelCase`.

Functions and class methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Names of functions and class methods use the same rules are variable names, with some
additions.

- If a method is used to get a property of an object that cannot be reset (example:
  `IceGrid::Mx()`), omit `get_` from the name.
- If a getter method has a corresponding setter method, their names should be
  *predictable*: `Foo::get_bar()` and `Foo::set_bar()`. In this case, *do not* omit `get_`
  from the name of the getter.

Namespaces
----------

Everything in PISM goes into the `pism` namespace. See the source code browser for more
namespaces (roughly one per sub-system).

Notable namespaces include:

- ``atmosphere``
- ``bed``
- ``calving``
- ``energy``
- ``frontalmelt``
- ``hydrology``
- ``ocean``
- ``rheology``
- ``sea_level``
- ``stressbalance``
- ``surface``

Using directives and declarations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Do *not* import all names from a namespace with `using namespace foo;`

Do import *specific* names with `using ::foo::bar;` in `.cc` files.


Formatting
----------

PISM includes a `.clang-format` file that makes it easy to re-format source to make it
conform to these guidelines.

To re-format a file, commit it to the repository, then run

.. code-block:: none

    clang-format -i filename.cc

(Here `-i` tells clang-format to edit files "in place." Note that editing in place is
safe because you added it to the repository.)

Logical operators should be surrounded by spaces
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: c++

   // Good
   if (a == b) {
     action();
   }

   // Bad
   if (a==b) {
     action();
   }

Commas and semicolons should be followed by a space
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: c++

   // Good
   function(a, b, c);

   // Bad
   function(a,b,c);
   function(a,b ,c);

Binary arithmetic operators should be surrounded by spaces
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: c++

   // Good
   f = x + y / (z * w);

   // Bad
   f = x+y/(z*w);

Statements
^^^^^^^^^^

One statement per line.

.. code-block:: c++

   // Good
   x = 0;
   y = x + 1;

   // Bad
   x = 0; y = x + 1;

Code indentation
^^^^^^^^^^^^^^^^

- Use two spaces per indentation level.
- Do not use tabs.
- Opening braces go with the keyword ("One True Brace Style").

Examples:

.. code-block:: c++

   int function(int arg) {
     return 64;
   }

   for (...) {
     something();
   }

   class Object {
   public:
     Object();
   };

Return statements
^^^^^^^^^^^^^^^^^

Return statements should appear on a line of their own.

Do not surround the return value with parenthesis if you don't have to.

.. code-block:: c++

   // Good
   int function(int argument) {
     if (argument != 0) {
       return 64;
     }
   }

   // Bad
   int function(int argument) {
     if (argument != 0) return(64);
   }


Conditionals
^^^^^^^^^^^^

- one space between `if` and the opening parenthesis
- no spaces between `(` and the condition (`(condition)`, not `( condition )`)
- all `if` blocks should use braces (`{` and `}`) *unless* it makes the code significantly
  harder to read
- `if (condition)` should always be on its own line
- the `else` should be on the same line as the closing parenthesis: `} else { ...`

.. code-block:: c++

   // Good
   if (condition) {
     action();
   }

   // Bad
   if (condition) action();

   // Sometimes acceptable:
   if (condition)
     action();

Loops
^^^^^

`for`, `while`, `do {...} unless` loops are formatted similarly to conditional blocks.

.. code-block:: c++

   for (int k = 0; k < N; ++k) {
     action(k);
   }

   while (condition) {
     action();
   }

   do {
     action();
   } unless (condition);

Error handling
--------------

First of all, PISM is written in C++, so unless we use a non-throwing `new` and completely
avoid STL, exceptions are something we have to live with. This means that we more or less
have to use exceptions to handle errors in PISM. (Mixing two error handling styles is a
*bad* idea.)

So, throw an exception to signal an error; PISM has a generic runtime error exception
class `pism::RuntimeError`.

To throw an exception with an informative message, use

.. code-block:: c++

   throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                 "format string %s", "data");

Error handling in a parallel program is hard. If all ranks in a communicator throw an
exception, that's fine. If some do and some don't PISM will hang as soon as one rank
performs a locking MPI call. I don't think we can prevent this in general, but we can
handle some cases.

Use

.. code-block:: c++

   ParallelSection rank0(communicator);
   try {
     if (rank == 0) {
       // something that may throw
     }
   } catch (...) {
     rank0.failed();
   }
   rank0.check();

to wrap code that is likely to fail on some (but not all) ranks. `rank0.failed()` prints
an error message from the rank that failed and `rank0.check()` calls `MPI_Allreduce(...)`
to tell other ranks in a communicator that everybody needs to throw an exception.
(`pism::ParallelSection::failed()` should be called in a `catch (...) {...}` block
**only**.)

In general one should not use `catch (...)`. It *should* be used in these three cases,
though:

1. With `pism::ParallelSection` (see above).
2. In callback functions passed to C libraries. (A callback is not allowed to throw, so we
   have to catch everything.)
3. In `main()` to catch all exceptions before terminating.

Performing an action every time a PISM exception is thrown
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The class `pism::RuntimeError` allows setting a "hook" that is called by the constructor
of `RuntimeError`. See the example below for a way to use it.

.. code-block:: c++

   #include <cstdio>

   #include "error_handling.hh"

   void hook(pism::RuntimeError *exception) {
     printf("throwing exception \"%s\"\n", exception->what());
   }

   int main(int argc, char **argv) {

     MPI_Init(&argc, &argv);

     pism::RuntimeError::set_hook(hook);

     try {
       throw pism::RuntimeError("oh no!");
     } catch (pism::RuntimeError &e) {
       printf("caught an exception \"%s\"!\n", e.what());
     }

     MPI_Finalize();

     return 0;
   }

Calling C code
^^^^^^^^^^^^^^

Check the return code of every C call and convert it to an exception if needed. Use macros
`PISM_CHK` and `PISM_C_CHK` for this.

When calling several C function in sequence, it may make sense to wrap them in a function.
Then we can check its return value and throw an exception if necessary.

.. code-block:: c++

   int call_petsc() {
     // Multiple PETSc calls here, followed by CHKERRQ(ierr).
     // This way we need to convert *one* return code into an exception, not many.
     return 0;
   }

   // elsewhere:
   int err = call_petsc(); PISM_C_CHK(err, 0, "call_petsc");

Use assert() for sanity-checks that should not be used in optimized code
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The `assert` macro should be used to check pre-conditions and post-conditions that can
fail *due to programming errors*.

**Do not** use `assert` to validate user input.

Note that *user input includes function arguments* for all functions and public members of
classes accessible using Python wrappers. (Use exceptions instead.)

Function design
---------------

Functions are the way to *manage complexity*. They are not just for code reuse: the main
benefit of creating a function is that a self-contained piece of code is easier both to
**understand** and **test**.

Functions should not have side effects (if at all possible). In particular, do not use and
especially do not *modify* "global" objects. If a function needs to modify an existing
field "in place", pass a reference to that field as an argument and document this argument
as an "input" or an "input/output" argument.

Function arguments; constness
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Pass C++ class instances by const reference *unless* an instance is modified in place.
This makes it easier to recognize *input* (read-only) and *output* arguments.

Do **not** use `const` when passing C types: `f()` and `g()` below are equivalent.

.. code-block:: c++

   double f(const double x) {
     return x*x;
   }

   double g(double x) {
     return x*x;
   }

Method versus a function
^^^^^^^^^^^^^^^^^^^^^^^^

Do **not** implement something as a class method if the same functionality can be
implemented as a standalone function. Turn a class method into a standalone function if
you notice that it uses the *public* class interface only.

Virtual methods
^^^^^^^^^^^^^^^

- Do not make a method virtual unless you have to.
- Public methods should not be virtual (create "non-virtual interfaces")
- **Never** add `__attribute__((noreturn))` to a virtual class method.

private versus protected versus public
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Most data members and class methods should be `private`.

Make it `protected` if it should be accessible from a derived class.

Make it `public` only if it is a part of the interface.
