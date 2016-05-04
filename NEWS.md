v0.8.3 (2015-05-04)
===================

- Remove BUNIT from reserved header keywords

v0.8.2 (2015-02-29)
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
