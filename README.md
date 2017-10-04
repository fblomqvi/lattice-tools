lattice-tools is a free collection of programs designed to make working with
lattices simple and efficient. However, at this moment, the tools are not ready
for use.

INSTALLATION
------------

In short:

     make
     make install


USAGE
-----

This package contains six different programs:

* __an-solver__ - solves the CVP on the lattice family A_n
* __lat-gen__ - generates random integer lattices
* __lat-sim__ - decoding simulations
* __lat-solve__ - solves the CVP on general lattices
* __lat-svp__ - solves the SVP for general lattices
* __rnd-point__ generates random points in Euclidian space

__lat-solve__ and __an-solve__ read the points to decode in binary format from
stdin and they are designed to be used together with __rnd-point__; the
intended usage is to pipe the output of __rnd-point__ to __lat-solve__ or
__an-solve__. __rnd-point__ can also be used as a stand-alone program.

__lat-gen__, __lat-sim__ and __lat-svp__ are intended as stand-alone programs.

See tha manual pages for these programs for more information.


NOTES
-----

__lat-svp__ only exists because future versions of this package will
need facilities for solving the SVP for internal use, and it would be stupid
not to spend 10 minutes on creating a program that exposes this functionality.
However, it only uses double precision arithmetic and will not be accurate in
all situations. Please use [fpLLL](http://github.com/fplll/fplll) or
[Magma](http://magma.maths.usyd.edu.au/magma/) if this is a concern.
