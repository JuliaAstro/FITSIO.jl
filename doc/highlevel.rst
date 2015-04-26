========================
High-level API reference
========================

This is currently an incomplete reference of methods for ``FITS``,
``FITSHeader``, and ``ImageHDU`` types.

File operations
---------------

.. function:: FITS(filename::String, mode::String="r")

   Open or create a FITS file. ``mode`` can be one of ``"r"``
   (read-only), ``"r+"`` (read-write) or ``"w"`` (write). In "write"
   mode, any existing file of the same name is overwritten.

.. function:: length(f::FITS)

   Return the number of HDUs in ``f``.

.. function:: getindex(f::FITS, i::Integer)

   Return the ``i``-th HDU. Same as ``f[i]``.

.. function:: getindex(f::FITS, name::String, ver::Int=0)

   Returns the HDU containing the given HDUNAME (EXTNAME) keyword,
   and optionally the given HDUVER (EXTVER) keyword.
   Same as ``f[name]`` or ``f[name, ver]``.

Header operations
-----------------

.. function:: readheader(hdu::HDU)

   Read the entire header from the given HDU and return a
   ``FITSHeader`` object. The value of each header record is parsed as
   ``Int``, ``Float64``, ``ASCIIString``, ``Bool`` or ``nothing``
   according to the FITS standard.  (If the value cannot be parsed
   according to the FITS standard, the value is stored as the raw unparsed
   ``ASCIIString``.)

Image operations
----------------

.. function:: read(hdu::ImageHDU)

   Read the entire image from disk.

.. function:: copy_section(hdu::ImageHDU, destination::FITS, r::Range...)

   Copy a rectangular section of an image and write it to a new FITS
   primary image or image extension. The new image HDU is appended to
   the end of the destination file; all the keywords in the input image
   will be copied to the output image. The common WCS keywords will be
   updated if necessary to correspond to the coordinates of the
   section. Examples:

   Copy the lower-left 200 x 200 pixel section of the image in ``hdu``
   to an open file, ``f``::
 
       copy_section(hdu, f, 1:200, 1:200)

   Same as above but only copy odd columns in y::

       copy_section(hdu, f, 1:200, 1:2:200)


Table Operations
----------------

.. function:: write(f::FITS, data::Dict)

   Create a new table extension and write data to it. If the FITS file is
   currently empty then a dummy primary array will be created before
   appending the table extension to it. ``data`` should be a dictionary
   with ASCIIString keys (giving the column names) and Array values
   (giving data to write to each column).

.. function:: read(hdu::TableHDU, colname)

   Read a column as an array from the given table HDU.

   The column name may contain wild card characters (*, ?, or #). The
   `*' wild card character matches any sequence of characters
   (including zero characters) and the `?' character matches any
   single character. The # wildcard will match any consecutive string
   of decimal digits (0-9). The string must match a unique column.
