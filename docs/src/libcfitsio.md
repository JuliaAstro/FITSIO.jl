# Libcfitsio submodule

```@meta
CurrentModule = FITSIO.Libcfitsio
```

The `Libcfitsio` submodule provides an interface familiar to users of the
[CFITSIO](http://heasarc.gsfc.nasa.gov/fitsio/) C library. It can be used with

```julia
using FITSIO.Libcfitsio
```

The functions exported by this module operate on `FITSFile` objects,
which is a thin wrapper around a pointer to a CFITSIO `fitsfile`.  For
the most part, the functions are thin wrappers around the CFITSIO
routines of the same names. Typically, they:

* Convert from Julia types to C types as necessary.
* Check the returned status value and raise an appropriate exception if
  non-zero.

!!! warning

    Note that these functions do not check if the file is still open
    before trying to access it. A segmentation fault can result from
    trying to operate on a closed file. (The main FITSIO interface
    always checks if the file is open before any operation.)


## File access

```@docs
fits_create_file
fits_clobber_file
fits_open_file
fits_open_table
fits_open_image
fits_open_data
fits_close_file
fits_delete_file
fits_file_name
```


## HDU Routines

The functions described in this section change the current
HDU and to find their number and type. The following is a short
example which shows how to use them:

```julia
num = fits_get_num_hdus(f)
println("Number of HDUs in the file: ", num)

for i = 1:num
    hdu_type = fits_movabs_hdu(f, i)
    println(i, ") hdu_type = ", hdu_type)
end
```

```@docs
fits_get_num_hdus
fits_movabs_hdu
fits_movrel_hdu
fits_movnam_hdu
```


## Header Keyword Routines

```@docs
fits_get_hdrspace
fits_read_keyword
fits_read_record
fits_read_keyn
fits_write_key
fits_write_record
fits_delete_record
fits_delete_key
fits_hdr2str
```


## Image HDU Routines

```@docs
fits_get_img_size
fits_create_img
fits_write_pix
fits_read_pix
```


## Table Routines

There are two functions to create a new HDU table extension:
`fits_create_ascii_table` and `fits_create_binary_table`. In general,
one should pick the second as binary tables require less space on the
disk and are more efficient to read and write. (Moreover, a few
datatypes are not supported in ASCII tables). In order to create a
table, the programmer must specify the characteristics of each column
by passing an array of tuples. Here is an example:

```julia
f = fits_create_file("!new.fits")
coldefs = [("SPEED", "1D", "m/s"),
           ("MASS", "1E", "kg"),
           ("PARTICLE", "20A", "Name")]
fits_create_binary_tbl(f, 10, coldefs, "PARTICLE")
```  

This example creates a table with room for 10 entries, each of them
describing the characteristics of a particle: its speed, its mass, and
its name (codified as a 20-character string). See the documentation of
`fits_create_ascii_tbl` for more details.

```@docs
fits_create_ascii_tbl
fits_create_binary_tbl
fits_get_coltype
fits_insert_rows
fits_delete_rows
fits_read_col
fits_write_col
```
