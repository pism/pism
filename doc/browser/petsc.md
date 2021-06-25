PETSc: An overview for PISM users {#petscID}
=====

## Preamble

In v0.5, because of abstractions especially in IceModelVec and its derived
classes, programmers building derived classes of IceModel now need minimal
knowledge of PETSc and especially PETSc DMs.  What *is* important to know is
that solving the SSA in parallel is relatively deep technology.  Very
specifically, because effective parallel iterative linear algebra requires
preconditioning, *the solution of the SSA stress balance depends on the number
of processors*.  This will also be true of high resolution parallel solution
of all other membrane-including stress balances (Blatter, Stokes).

## Introduction

The PETSc library [@ref petsc-user-ref] provides essential support for distributed 
arrays and linear solvers in a parallel computing environment.  "PETSc" stands for 
"Portable, Extensible Toolkit for Scientific Computation."  It is a suite of data structures 
and routines in C for the scalable parallel solution of partial differential equations and their 
scientific applications.  Large parts of PETSc relate especially to finite-difference-type regular, 
rectangular grids but PETSc has been used for unstructured grids as well.

PETSc employs the MPI standard for all message-passing communication but it is deliberately at a 
higher level of abstraction than MPI.  We find that using PETSc protects the programmer from 
explicit consideration of message passing about 90% of the time.

Documentation for PETSc is at 

http://www.mcs.anl.gov/petsc/.

PISM is a C++ program and PETSc is a C library, but all PISM calls to PETSc use
the PETSc C API, not the C++ API.


### PETSc types

Most important variables deep inside PISM are of these PETSc types:
  - `DM`
  - `Vec`
  - `Mat`
  - `KSP`

In fact most of the PETSc types merely declare pointers, as these are
abstract data types.  Variables must be created with calls to functions like `DMDACreate2d()`, 
`VecCreate()`, etc., and destroyed when no longer needed.  These concepts are
largely buried inside IceGrid (esp. DM type) and IceModelVec (esp. Vec type).
But Mat and KSP types are exposed in the IceModel code which solves the SSA.


### Distributed arrays and vectors

PETSc has an abstract date type called a Distributed Array, which is a kind of
a Distributed Mesh (type `DM`). `DMs` contain information about the grid and
stencil. They set up *ghosted* values which are intended to duplicate the
values on neighboring processors so communication can be done in big batches.

Vectors (type `Vec`) are created *using* a `DM` with `DMCreateLocalVector()` and similar 
procedures.  These vectors are distributed across the processors according to the information
in the `DM`.  Just for emphasis: *PISM distributes itself across the processors by calling
PETSc procedures which create a `DM`*.

There are two parameters of note when creating a `DMDA`: *stencil type* and 
*stencil width.*  

The stencil types are `DMDA_STENCIL_STAR` and `DMDA_STENCIL_BOX`.  
They are generalizations of the five point and nine point stencils typical of two-dimensional 
discretizations of second order partial derivatives.  Using `DMDA_STENCIL_STAR` means
that ghosted points are needed only along the coordinate axes, while `DMDA_STENCIL_BOX` 
indicates that ghosted points are needed in the box-shaped %region surrounding each point, 
and thus each processor.  (On processor *N,* information owned by a neighboring processor 
but duplicated onto *N* is called *ghosted.*  A picture is worth a thousand words here, so 
find the appropriate picture in the PETSc manual!)  We only need `BOX` style stencils when gradient 
terms must be evaluated on a staggered grid (*h* in SIA regions and @f$\bar{u},\bar{v}@f$
in the computation of effective viscosity in SSA regions).  Keeping all other 
two-dimensional vectors on a `STAR` type stencil would reduce the necessary communication
slightly, but for now all two dimensional arrays are based on `BOX`.

The stencil width indicates how many points in each direction will be needed.  We never 
need a stencil width greater than 1.

The three dimensional distributed arrays are distributed across processors in the same way
as two-dimensional arrays, in the sense that *each column of the ice (or bedrock in that 
thermal model) is owned by exactly one processor*.  From the point of view of parallelization,
our problem is two-dimensional.  All three-dimensional arrays have stencil type `STAR` in the 
horizontal directions.

One point of confusion (we admit ...) is the redefinition of the @f$x,y,z@f$ axes. Contrary 
to the PETSc default, our @f$z@f$ axis changes most rapidly through memory while the
@f$x@f$ axis changes most slowly.  The desirable consequence, however, is that our C
arrays are addressed as `u[i][j][k]` where `i,j,k` are the coordinate 
indices in the directions @f$x,y,z@f$ respectively.


### Accessing the processor's part of a DM-managed Vec

`DM` -based vectors must be accessed by a call to `DMDAVecGetArray()` before the access and
a call to `DMDAVecRestoreArray()` after.  This just means that one gets a valid pointer 
to the actual memory, for the abstract data type `Vec`.  The point has type 
`double**` for a 2-dimensional array and type `double***` for a 3-dimensional 
array.

The resulting pointer can be addressed using normal (to the extent that C is "normal" in 
this regard ...) multidimensional array indexing.  Furthermore, the indices themselves 
allow one to pretend that the given processor is addressing the entire array, even though only a
chunk is stored on the local processor.  *No message passing occurs here.*  Instead, the crucial
idea is that the `DM` knows what are the valid index ranges for each processor.  That's why 
essentially all loops over two dimensional arrays in PISM look like this:


      for (PetscInt i = grid.xs; i < grid.xs + grid.xm; ++i) {
        for (PetscInt j = grid.ys; j < grid.ys + grid.ym; ++j) {
        ... vH(i,j) ...
        }
      }


The integers `grid.xs`, `grid.xm`, `grid.ys`, `grid.ym`, are just the numbers that come
from a call to `DMDAGetLocalInfo()`, just after the `DM` is created, but normally
a programmer just calls IceGrid.


### "Local" versus "global" Vecs

PETSc `DM` -based vectors can be "local" or "global".   (This weird PETSc language should
probably be translated as "has ghosts" or "sans ghosts".)  "Local" vectors include 
space for the ghosted points.  That is, when `DAVecGetArray()` is called, the resulting 
array can be indexed on all of the local portion owned by the processor, *plus* the 
ghosted points.  The processor owns *copies* of the values at grid points owned by the neighboring
processors, according to the stencil type and stencil width in the `DA`.  The processor should 
treat these ghost values as "read-only", because of what happens at the communication stage 
(which we describe in a moment).

*Example.*  We create a local `Vec` `v` which we want to have `STAR` type stencil.  Assuming 
`da` is `DM` created with `DMDA_STENCIL_STAR` then we do this to create `v` and get an array 
pointer `arr` into the local portion of it:

    DMDACreateLocalVector(da, &v);
    DMDAVecGetArray(da, v, &arr);

Now the processor can access `arr[i][j]` for every 

    grid.xs <= i <= grid.xs + grid.xm - 1

and

    grid.ys <= j <= grid.ys + grid.ym - 1

But this is just the regular (owned by the processor) portion.  In addition, all of these 
are valid for every `i` and `j` in the above ranges:

    arr[i+1][j],   arr[i-1][j],   arr[i][j+1],   arr[i][j-1].

(*If* the `da` were of type `DMDA_STENCIL_BOX` then, in addition 

    arr[i+1][j+1],   arr[i-1][j+1],   arr[i+1][j-1],   arr[i-1][j-1].

would all be valid.)  Once the processor is done with modifying its portion of `v` (i.e. not including 
the ghosts, although it may read them) then we would want to update the ghost values so that
all processors agree on the values at those grid points:

    DMDAVecRestoreArray(da, v, &arr);
    DMDALocalToLocalBegin(da, v, INSERT_VALUES, v);
    DMDALocalToLocalEnd(da, v, INSERT_VALUES, v);

(Note we release the array pointer before the communication stage, although that may not be
essential.)

`DMDALocalToLocalBegin()` and then `DMDALocalToLocalEnd()` are called to update (communicate) the 
ghosted values before the next stage at which they will be needed.  This pair of commands 
perform a stage of actual (MPI) message passing.  Only enough values are passed around to 
update the ghosts.  The size of the messages is relatively small, and *latency is more important
(in slowing performance) than bandwidth*.  The significance of `STAR` versus `BOX` is not
in the size of the extra memory, but in the extra messages which must be passed to the "diagonal"
neighboring processors.

Global vectors do not hold ghosted values.  Certain array operations are more appropriate to
these kind of vectors, like viewing arrays in graphics windows.  (Of course, viewing an array 
distributed across all the processors requires message passing.  I think PETSc sends all the
values to processor 0 to display then.)  In any case, when ghosts are not needed there is
no need to allocate space for them, and a "global" `Vec` is appropriate.  Local vectors
typically need to be mapped to global vectors before viewing or using in a linear system, but
this is an entirely "local" operation which does not require message passing (`DMLocalToGlobalBegin/End()`).

At risk of repetition, most of the above distinction between "local" and "global"
is buried in IceModelVec abstraction.

### Solving linear systems

PETSc is designed for solving large, sparse systems in parallel.  Iterative methods are 
the name of the game and especially Krylov subspace methods such as conjugate gradients and 
GMRES [@ref TrefethenBau, @ref Saad].  For consistency, all methods use the same 
Krylov-subspace-method-with-preconditioner `KSP` interface, even though some are really direct 
methods.

The place in PISM where this is already used is in the solution of the linearized SSA velocity
problem.  See the documentation for SSAFD and SSAFEM in this *Manual.*

One starts by declaring an object of type `KSP`; in PISM this is done in IceModel::createVecs().
Various options can be set by the program, but, if the program calls `KSPSetFromOptions()` 
before any linear systems are solved, which PISM does, then user options like 
`-ksp_type`, `-pc_type`, and `-ksp_rtol` can be read to control (respectively) 
which Krylov method is used, which preconditioner is used, and what relative 
tolerance will be used as the convergence criterion of the Krylov solver.  The default is
GMRES(30) with ILU preconditioning.

(As a general mechanism, PETSc has an options database which holds command line options.
All PISM options use this database.  See IceModel::setFromOptions().)

To solve the system, a matrix must be attached to the `KSP`.  The first time
`KSPSolve()` is called, the matrix will be factored by the preconditioner and reused
when the system is called for additional right hand sides.  In the case of PISM and the 
solution of the balance of momentum equations for the SSA, the matrix changes because the
equations are nonlinear (because the effective viscosity depends on the strain rates).

The default matrix format is similar to the *Matlab* `sparse` format.  Each processor
owns a range of rows.  Elements in matrices and vectors can be set using `MatSetValues()`
and `VecSetValues()`.  These routines use a "global" indexing which is not based on the regular
grid and stencil choices in the `DA`.  One can set any value on any processor but a lot of
message passing has to occur; fortunately PETSc manages that all.  Values are
cached until one calls `MatAssemblyBegin()` followed by `MatAssemblyEnd()` to
communicate the values.

### Additional PETSc utility functions

Quite ellaborate error tracing and performance monitoring is possible with PETSc.  All
functions return `PetscErrorCode` which should be checked by the macro `CHKERRQ()`.
Therefore runtime errors normally print sufficient traceback information.
If this information is not present (because the error did not get traced back), you 
may need to use a debugger which is accessible with the command line (PETSc) options 
`-start_in_debugger` and `-on_error_attach_debugger`.  Option `-log_summary`
is useful to get some diagnostics written to the terminal.

The `PetscViewer` interface allows PETSc objects to be displayed. One
can "view" a `Vec` to a graphical (X) window, but the "view" can be
saving the `Vec` to a binary file on disk, in plain text to the
terminal (standard out), or even by sending to a running instance of
*Matlab.* In PISM the X views are available using the `-view`
option. See the "User's Manual". When viewing `Vec` s graphically in
multiprocessor jobs, the display may have to be set on the command
line, for instance as `-display :0` or similar; this must be given as
the final option. For example,

    mpiexec -n 2 pismr -eisII A -view thk -display :0

allows a two processor run to view the ice thickness.
