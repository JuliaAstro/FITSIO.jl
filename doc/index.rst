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

