=======================
FITSIO.jl Documentation
=======================

A Julia_ module for reading and writing Flexible Image Transport
System (FITS) files, based on the CFITSIO_ library.

The high-level interface is inspired by Erin Sheldon's FITSIO_ python module.

.. _Julia: http://julialang.org
.. _CFITSIO: http://heasarc.gsfc.nasa.gov/fitsio/
.. _FITSIO: https://github.com/esheldon/fitsio

Install
-------

::

    julia> Pkg.add("FITSIO")

The cfitsio library is automatically downloaded and compiled. This
will not interfere with any other versions of cfitsio on your system.

Documentation
-------------

.. toctree::
   :maxdepth: 1

   usage
   highlevel
   lowlevel
