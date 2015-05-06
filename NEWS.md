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
