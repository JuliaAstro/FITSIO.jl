========================
High-level API reference
========================

This is currently an incomplete reference of methods for ``FITS``,
``FITSHeader``, and ``ImageHDU`` types.

FITS operations
---------------

.. function:: FITS(filename::String, mode::String="r")

   Open or create a FITS file. ``mode`` can be one of ``"r"``
   (read-only), ``"r+"`` (read-write) or ``"w"`` (write). In "write"
   mode, any existing file of the same name is overwritten.

.. function:: length(f::FITS)

   Return the number of HDUs in ``f``.

.. function:: getindex(f::FITS, i::Integer)
.. function:: getindex(f::FITS, name::String, ver::Int=0)

   In the first form, returns the ``i``-th HDU. Same as ``f[i]``.
   In the second form returns the HDU by HDUNAME (EXTNAME), and optionally
   HDUVER (EXTVER). Same as ``f[name]`` or ``f[name, ver]``.

Header operations
-----------------

.. function:: readheader(hdu::HDU)

   Read the entire header from the given HDU and return a
   ``FITSHeader`` object.

Image operations
----------------

.. function:: read(hdu::ImageHDU)

   Read the entire image from disk.

