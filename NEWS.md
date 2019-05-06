v0.14.0 (2019-05-06)
====================

- Julia v1.0 is required.

v0.13.0 (2019-01-11)
====================

- Upgrade CFITSIO library to v3.45.

v0.12.0 (2018-09-18)
====================

- Drop support for Julia 0.5.
- Column names in FITS files are now case sensitive by default.  Functions
  `read(hdu, colname)` and `Libcfitsio.fits_get_colnum` have a new optional
  keyword, `case_sensitive`, which is `true` by default.

v0.11.0 (2017-09-28)
====================

- Drop support for Julia 0.4.
- Methods now accept `String` and, for strings written to the FITS file,
  an error is raised for non-ASCII strings.
- New `write_key()` function for writing to an HDU header.
- Homebrew.jl used for libcfitsio dependency on OS X.
- Remove extra newline from `show` output.
- Fix some deprecation warnings on Julia 0.7.
- Inline docstrings and HTML documentation generated with Documenter.jl.

v0.10.0 (2017-06-13)
====================

- Add `FITS(::Vector{UInt8})` constructor: support for memory-backed files.
- Add support for `FITS("file.fits") do f; ...; end`, which automatically
  closes the file at the end of the block.
- Add Libcfitsio functions: `fits_read_key_str`, `fits_read_key_lng`,
  `fits_read_keys_lng`, `fits_write_date` low-level functions.
- Error message from opening a file now includes file name.

v0.9.0 (2017-01-24)
===================

- Fix deprecation warnings on Julia 0.6.
- Drop support for Julia 0.3.

v0.8.4 (2016-08-09)
===================

- Fix deprecation warnings on v0.5

v0.8.3 (2016-05-04)
===================

- Remove BUNIT from reserved header keywords

v0.8.2 (2016-02-29)
===================

- enable precompilation

v0.8.1 (2015-09-30)
===================

- Fix deprecation warnings on v0.4

v0.8.0 (2015-07-03)
===================

## New features

- Windows support (32-bit & 64-bit), thanks to help from Tony Kelman
  (requires master of BinDeps)
- New method `read_header(hdu, ASCIIString)` returns entire header as
  a single string.
- bump cfitsio version from 3.360 to 3.370

## Bug fixes

- fix `show()` for an empty FITS file.
- fix several issues with `Clong` on 64-bit Windows.

v0.7.0 (2015-05-31)
===================

## New Features & Improvements

- Read and write variable length columns in binary tables.
- Iterate over HDUs in a file with `for hdu in f; ...; end`
- Improved `show()` methods give more and better information
  for `FITS` and `HDU` types.
- Write name and version of HDU when creating a new HDU with `name` and
  `ver` keywords.
- New function `FITSIO.libcfitsio_version()` (unexported) returns CFITSIO
  library version.

## Deprecations

- `Libcfitsio.fits_read_num_rowsll` deprecated.
  Use `Libcfitsio.fits_read_num_rows`.


v0.6.0 (2015-05-06)
===================

## New Features

- Read and write table extensions (both ASCII and binary) in
  high-level API.

## Breaking changes

The low-level API functions (starting with `fits_*`) have been moved to
the `Libcfitsio` sumodule. If you are using these functions, simply add

```julia
using FITSIO.Libcfitsio
```

in place of, or in addition to, `using FITSIO`.

## Deprecations

- `readkey` renamed to `read_key`
- `readheader` renamed to `read_header`
- `getcomment` renamed to `get_comment`
- `setcomment!` renamed to `set_comment!`
- `hdu[i:j, :]` replaced by `read(hdu, i:j, :)` for reading subsets
  of image extensions.

In the low-level interface (now in `Libcfitsio`):

- `fits_get_col_repeat` deprecated. Use `fits_get_coltype`, which
  returns the column typecode in addition to width and repeat values.

- Methods of `fits_read_col` and `fits_write_col` specifying the type
  explicitly are deprecated because the type is redundant. Simply
  remove the explicit type argument.

## Internal

- Cleanup and correction of type specifications (e.g., `Cint` in place
  of `Int32`, `ASCIIString` in place of `String`)
