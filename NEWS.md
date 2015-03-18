v0.6.0 (unreleased)
===================

## Deprecations in low-level c-style interface

- `fits_get_col_repeat` deprecated. Use `fits_get_coltype`, which
  returns the column typecode in addition to width and repeat values.

- Methods of `fits_read_col` and `fits_write_col` specifying the type
  explicitly are deprecated because the type is redundant. Simply
  remove the explicit type argument.

## Internal

- Cleanup and correction of type specifications (e.g., `Cint` in place
  of `Int32`, `ASCIIString` in place of `String`)