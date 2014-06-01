FITSIO.jl
=========

Flexible Image Transport System (FITS) support for Julia

Installation
------------

```jlcon
julia> Pkg.add("FITSIO")
```

Simple Example
--------------

```jlcon

julia> using FITSIO

julia> f = FITS("file.fits", "r")
file: test.fits
mode: r
extnum exttype         extname
1      image_hdu
```

Complete Documentation
----------------------

https://julia-fitsio.readthedocs.org/

