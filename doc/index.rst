FITSIO --- FITS File I/O
========================

.. module:: FITSIO
   :synopsis: Read and write FITS files.

A wrapper for the CFITSIO_ library. 

.. _CFITSIO: http://heasarc.gsfc.nasa.gov/fitsio/

File Access Routines
--------------------

.. function:: fits_create_file(filename::String)

   Create and open a new empty output FITS file.

.. function:: fits_clobber_file(filename::String)

   Like fits_create_file, but overwrites ``filename`` if it exists.

.. function:: fits_open_file(filename::String)

   Open an existing data file.

.. function:: fits_open_table(filename::String)

   Open an existing data file (like :func:`fits_open_file`) and move
   to the first HDU containing either an ASCII or a binary table.

.. function:: fits_open_image(filename::String)

   Open an existing data file (like :func:`fits_open_file`) and move
   to the first HDU containing an image.

.. function:: fits_open_data(filename::String)

   Open an existing data file (like :func:`fits_open_file`) and move
   to the first HDU containing either an image or a table.

.. function:: fits_close_file(f::FITSFile)

   Close a previously opened FITS file.

.. function:: fits_delete_file(f::FITSFile)

   Close an opened FITS file (like :func:`fits_close_file`) and
   removes it from the disk

.. function:: fits_file_name(f::FITSFile)

   Return the name of the file associated with object `f`.

HDU Routines
------------

The functions described in this section allow to change the default
HDU and to find their number and type. The following is a short
example which shows how to use them:

.. code-block:: julia

  num = fits_get_num_hdus(f)
  println("Number of HDUs in the file: ", num)

  for i = 1:num
      hdu_type = fits_movabs_hdu(f, i)
      println(i, ") hdu_type = ", hdu_type)
  end


.. function:: fits_get_num_hdus(f::FITSFile)

   Return the number of HDUs in the file.

.. function:: fits_movabs_hdu(f::FITSFile, hduNum::Integer)

   Change the current HDU to the value specified by `hduNum`, and
   return a symbol describing the type of the HDU. Possible symbols
   are: ``:image_hdu``, ``:ascii_table``, or ``:binary_table``.

   The value of `hduNum` must range between 1 and
   the value returned by :func:`fits_get_num_hdus`.

.. function:: fits_movrel_hdu(f::FITSFile, hduNum::Integer)

   Change the current HDU by moving forward or backward by `hduNum`
   HDUs (positive means forward), and return the same as
   :func:`fits_movabs_hdu`.


Header Keyword Routines
-----------------------

.. function:: fits_get_hdrspace(f::FITSFile) -> (keysexist, morekeys)

   Return the number of existing keywords (not counting the END keyword)
   and the amount of space currently available for more keywords.

.. function:: fits_read_keyword(f::FITSFile, keyname::String) -> (value, comment)

   Return the specified keyword.

.. function:: fits_read_record(f::FITSFile, keynum::Int) -> String

   Return the nth header record in the CHU. The first keyword in the header is at ``keynum = 1``.

.. function:: fits_read_keyn(f::FITSFile, keynum::Int) -> (name, value, comment)

   Return the nth header record in the CHU. The first keyword in the header is at ``keynum = 1``.

.. function:: fits_write_key(f::FITSFile, keyname::String, value, comment::String)

   Write a keyword of the appropriate data type into the CHU.

.. function:: fits_write_record(f::FITSFile, card::String)

   Write a user specified keyword record into the CHU.

.. function:: fits_delete_record(f::FITSFile, keynum::Int)

   Delete the keyword record at the specified index.

.. function:: fits_delete_key(f::FITSFile, keyname::String)

   Delete the keyword named ``keyname``.

Primary Array Routines
----------------------

.. function:: fits_get_img_size(f::FITSFile)

   Get the dimensions of the image.

.. function:: fits_create_img(f::FITSFile, t::Type, naxes::Vector{Int})

   Create a new primary array or IMAGE extension with a specified data type and size.

.. function:: fits_write_pix(f::FITSFile, fpixel::Vector{Int}, nelements::Int, data::Array)

   Write pixels from `data` into the FITS file.

.. function:: fits_read_pix(f::FITSFile, fpixel::Vector{Int}, nelements::Int, data::Array)

   Read pixels from the FITS file into ``data``.


Table Routines
--------------

To create ASCII/binary tables in a new HDU, the FITSIO.jl library
provides two functions: :func:`fits_create_ascii_table` and
:func:`fits_create_binary_table`. In general, one should pick the
second as binary tables require less space on the disk and are more
efficient to read and write. (Moreover, a few datatypes are not
supported in ASCII tables). In order to create a table, the programmer
must specify the characteristics of each column by passing an array of
tuples. See the documentation of :func:`fits_create_ascii_table` for
more details.

Here is an example:

.. code-block:: julia

   f = fits_create_file("!new.fits")
   coldefs = [("SPEED", "1D", "m/s"),
              ("MASS", "1E", "kg"),
              ("PARTICLE", "20A", "Name")]
   fits_create_binary_tbl(f, 10, coldefs, "PARTICLE")
  

This example creates a table with room for 10 entries, each of them
describing the characteristics of a particle: its speed, its mass, and
its name (codified as a 20-character string).

.. function:: fits_create_ascii_table(f::FITSFile, numrows::Integer, coldefs::Array{ColumnDef}, extname::String)

   Append a new HDU containing an ASCII table. The table will have
   `numrows` rows (this parameter can be set to zero), each
   initialized with the default value. The columns are specified by
   the `coldefs` variable, which is an array of tuples. Each tuple
   must have three string fields:

   1. The name of the column.
   2. The data type and the repetition count. It must be a string made
      by a number (the repetition count) followed by a letter
      specifying the type (in the example above, ``D`` stands for
      `Float64`, ``E`` stands for Float32, ``A`` stands for ``Char``).
      Refer to the CFITSIO documentation for more information about
      the syntax of this parameter.
   3. The measure unit of this field. This is used only as a comment.

   The value of `extname` sets the "extended name" of the
   table, i.e., a string that in some situations can be used to refer
   to the HDU itself.

   Note that, unlike for binary tables, CFITSIO puts some limitations
   to the types that can be used in an ASCII table column. Refer to
   the CFITSIO manual for further information.

   See also :func:`fits_create_binary_table` for a similar function
   which creates binary tables.

.. function:: fits_create_binary_table(f::FITSFile, numrows::Integer, coldefs::Array{ColumnDef}, extname::String)

   Append a new HDU containing a binary table. The meaning of the
   parameters is the same as in a call to
   :func:`fits_create_ascii_table`.

.. function:: fits_get_col_repeat(f::FITSFile, colnum::Integer)

   Provided that the current HDU contains either an ASCII or binary
   table, this function returns a tuple containing two elements:

   1. the repetition count for the column at position `colnum`
      (starting from 1), and
   2. the optimal number of characters needed to print the value of
      any field contained in this column.

.. function:: fits_insert_rows(f::FITSFile, firstrow::integer, nrows::Integer)

   Insert a number of rows equal to `nrows` after the row number
   `firstrow`. The elements in each row are initialized to their
   default value: you can modify them later using
   :func:`fits_write_col`.

   Since the first row is at position 1, in order to insert rows
   *before* the first one `firstrow` must be equal to zero.

   See also :func:`fits_delete_rows`.

.. function:: fits_delete_rows(f::FITSFile, firstrow::integer, nrows::Integer)

   Delete `nrows` rows, starting from the one at position `firstrow`
   (the first row has index 1).

   See also :func:`fits_insert_rows`.


.. function:: fits_read_col{T}(f::FITSFile, ::Type{T}, colnum::Int, firstrow::Int64, firstelem::Int64, data::Array{T})

   Read data from one column of an ASCII/binary table and convert the
   data into the specified type `T`. The column number is specified by
   *colnum* (the first column has ``colnum=1``). The elements to be
   read start from the row number `firstrow`; in case each cell
   contains more than one element (i.e., the "repetition count" of the
   field is greater than one), `firstelem` allows to specify which is
   the first element to be read. The overall number of elements is
   specified by the length of the array `data`, which at the end of
   the call will be filled with the elements read from the column.

.. function:: fits_write_col{T}(f::FITSFile, ::Type{T}, colnum::Int, firstrow::Int64, firstelem::Int64, data::Array{T})

   Write some data in one column of a ASCII/binary table. The column
   number is specified by *colnum* (the first column has
   ``colnum=1``). The first element is written at the position
   `firstelem` within the row number `firstrow` (both the indexes
   start from one).

   If there is no room for the elements, new rows will be created. (It
   is therefore useless to call :func:`fits_insert_rows` if you only
   need to *append* elements to the end of a table.)
