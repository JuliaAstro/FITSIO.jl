FITSIO.jl: Flexible Image Transport System (FITS) support for Julia
===================================================================

Installation
------------

```jlcon
julia> Pkg.add("FITSIO")
```

Requires the [cfitsio](http://heasarc.gsfc.nasa.gov/fitsio/) library
to be installed on your system.

Simple Example
--------------

```jlcon
julia> require("FITSIO")
julia> using FITSIO
julia> x = fitsread("/path/to/file.fits")
```

Documentation
-------------

* [Low-level API](https://julia-fitsio.readthedocs.org/)

