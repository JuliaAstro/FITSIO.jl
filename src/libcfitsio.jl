# cfitsio.jl - C-style interface to CFITSIO functions:
#
# - Function names closely mirror the C interface (e.g., `fits_open_file()`).
# - Functions operate on `FITSFile`, a thin wrapper for `fitsfile` C struct
#   (`FITSFile` has concept of "current HDU", as in CFITSIO).
# - Note that the wrapper functions *do* check the return status from CFITSIO
#   and throw an error with the appropriate message.
#
#
# The following table gives the correspondances between CFITSIO "types",
# the BITPIX keyword and Julia types.
#
#     -------------------------------------------------
#     CODE  CFISTIO         Julia     Comments
#     -------------------------------------------------
#           int             Cint
#           long            Clong
#           LONGLONG        Int64     64-bit integer
#     -------------------------------------------------
#     -------- FITS BITPIX ----------------------------
#        8  BYTE_IMG        Uint8
#       16  SHORT_IMG       Int16
#       32  LONG_IMG        Int32
#       64  LONGLONG_IMG    Int64
#      -32  FLOAT_IMG       Float32
#      -64  DOUBLE_IMG      Float64
#      -------- cfitsio "aliases" ---------------------
#       10  SBYTE_IMG       Int8     written as: BITPIX = 8, BSCALE = 1,
#                                                BZERO = -128
#       20  USHORT_IMG      Uint16   written as: BITPIX = 16, BSCALE = 1,
#                                                BZERO = 32768
#       40  LONG_IMG        Uint32   written as: BITPIX = 32, BSCALE = 1,
#                                                BZERO = 2147483648
#     -------------------------------------------------
#     -------- FITS TABLE DATA TYPES ------------------
#        1  TBIT
#       11  TBYTE           Cuchar = Uint8
#       12  TSBYTE          Cchar = Int8
#       14  TLOGICAL        Bool
#       16  TSTRING         String
#       20  TUSHORT         Cushort
#       21  TSHORT          Cshort
#       30  TUINT           Cuint
#       31  TINT            Cint
#       40  TULONG          Culong
#       41  TLONG           Clong
#       42  TFLOAT          Cfloat
#       81  TLONGLONG       Int64
#       82  TDOUBLE         Cdouble
#       83  TCOMPLEX        Complex{Cfloat}
#      163  TDBLCOMPLEX     Complex{Cdouble}
#     -------------------------------------------------
#

__precompile__()


module Libcfitsio

export FITSFile,
       FITSMemoryHandle,
       fits_assert_open,
       fits_clobber_file,
       fits_close_file,
       fits_copy_image_section,
       fits_create_ascii_tbl,
       fits_create_binary_tbl,
       fits_create_file,
       fits_create_img,
       fits_delete_file,
       fits_delete_key,
       fits_delete_record,
       fits_delete_rows,
       fits_file_mode,
       fits_file_name,
       fits_get_hdrspace,
       fits_get_hdu_num,
       fits_get_hdu_type,
       fits_get_img_dim,
       fits_get_img_equivtype,
       fits_get_img_size,
       fits_get_img_type,
       fits_get_num_cols,
       fits_get_num_hdus,
       fits_get_num_rows,
       fits_get_rowsize,
       fits_get_colnum,
       fits_get_coltype,
       fits_get_eqcoltype,
       fits_get_version,
       fits_read_tdim,
       fits_hdr2str,
       fits_insert_rows,
       fits_movabs_hdu,
       fits_movrel_hdu,
       fits_movnam_hdu,
       fits_open_data,
       fits_open_file,
       fits_open_image,
       fits_open_table,
       fits_open_memfile,
       fits_read_col,
       fits_read_descript,
       fits_read_keyn,
       fits_read_key_str,
       fits_read_key_lng,
       fits_read_keys_lng,
       fits_read_keyword,
       fits_read_pix,
       fits_read_record,
       fits_read_subset,
       fits_update_key,
       fits_write_col,
       fits_write_date,
       fits_write_comment,
       fits_write_history,
       fits_write_key,
       fits_write_pix,
       fits_write_record,
       fits_write_tdim

if isfile(joinpath(dirname(@__FILE__),"..","deps","deps.jl"))
    include("../deps/deps.jl")
else
    error("FITSIO not properly installed. Please run Pkg.build(\"FITSIO\")")
end

const TYPE_FROM_BITPIX = Dict{Cint, DataType}()
for (T, code) in ((UInt8,     8), # BYTE_IMG
                  (Int16,    16), # SHORT_IMG
                  (Int32,    32), # LONG_IMG
                  (Int64,    64), # LONGLONG_IMG
                  (Float32, -32), # FLOAT_IMG
                  (Float64, -64), # DOUBLE_IMG
                  (Int8,     10), # SBYTE_IMG
                  (UInt16,   20), # USHORT_IMG
                  (UInt32,   40)) # ULONG_IMG
    local value = Cint(code)
    @eval begin
        TYPE_FROM_BITPIX[$value] = $T
        bitpix_from_type(::Type{$T}) = $value
    end
end

for (T, code) in ((UInt8,       11),
                  (Int8,        12),
                  (Bool,        14),
                  (String,      16),
                  (Cushort,     20),
                  (Cshort,      21),
                  (Cuint,       30),
                  (Cint,        31),
                  (Int64,       81),
                  (Float32,     42),
                  (Float64,     82),
                  (ComplexF32,  83),
                  (ComplexF64, 163))
    @eval cfitsio_typecode(::Type{$T}) = Cint($code)
end

# Above, we don't define a method for Clong because it is either Cint (Int32)
# or Int64 depending on the platform, and those methods are already defined.
# Culong is either UInt64 or Cuint depending on platform. Only define it if
# not already defined.
if Culong !== Cuint
    cfitsio_typecode(::Type{Culong}) = Cint(40)
end

# -----------------------------------------------------------------------------
# FITSFile type

mutable struct FITSFile
    ptr::Ptr{Cvoid}

    function FITSFile(ptr::Ptr{Cvoid})
        f = new(ptr)
        finalizer(f, fits_close_file)
        f
    end
end

# FITS wants to be able to update the ptr, so keep them
# in a mutable struct
mutable struct FITSMemoryHandle
    ptr::Ptr{Cvoid}
    size::Csize_t
end
FITSMemoryHandle() = FITSMemoryHandle(C_NULL, 0)

# -----------------------------------------------------------------------------
# error messaging

function fits_assert_open(f::FITSFile)
    if f.ptr == C_NULL
        error("attempt to access closed FITS file")
    end
end

function fits_get_errstatus(status::Cint)
    msg = Vector{UInt8}(undef, 31)
    ccall((:ffgerr,libcfitsio), Cvoid, (Cint,Ptr{UInt8}), status, msg)
    unsafe_string(pointer(msg))
end

function fits_assert_ok(status::Cint, filename = nothing)
    if status != 0
        prefix = filename == nothing ? "" : "While processing file `$filename`: "
        error(string(prefix, fits_get_errstatus(status)))
    end
end

fits_assert_isascii(str::String) =
    !isascii(str) && error("FITS file format accepts ASCII strings only")

fits_get_version() = ccall((:ffvers, libcfitsio), Cfloat, (Ref{Cfloat},), 0.)

# -----------------------------------------------------------------------------
# file access & info functions

"""
    fits_create_file(filename::AbstractString)

Create and open a new empty output `FITSFile`.
"""
function fits_create_file(filename::AbstractString)
    ptr = Ref{Ptr{Cvoid}}()
    status = Ref{Cint}(0)
    ccall((:ffinit,libcfitsio), Cint, (Ref{Ptr{Cvoid}},Ptr{UInt8},Ref{Cint}),
          ptr, filename, status)
    fits_assert_ok(status[], filename)
    FITSFile(ptr[])
end

"""
    fits_clobber_file(filename::AbstractString)

Like [`fits_create_file`](@ref), but overwrites `filename` if it exists.
"""
fits_clobber_file(filename::AbstractString) = fits_create_file("!"*filename)

"""
    fits_open_data(filename::String)

Open an existing data file (like [`fits_open_file`](@ref)) and move to the first HDU
containing either an image or a table.
"""
function fits_open_data end

"""
    fits_open_file(filename::String)

Open an existing data file.
"""
function fits_open_file end

"""
    fits_open_image(filename::String)

Open an existing data file (like [`fits_open_file`](@ref)) and move to the first
HDU containing an image.
"""
function fits_open_image end

"""
    fits_open_table(filename::String)

Open an existing data file (like [`fits_open_file`](@ref)) and move to the first
HDU containing either an ASCII or a binary table.
"""
function fits_open_table end

for (a,b) in ((:fits_open_data, "ffdopn"),
              (:fits_open_file, "ffopen"),
              (:fits_open_image,"ffiopn"),
              (:fits_open_table,"fftopn"))
    @eval begin
        function ($a)(filename::AbstractString, mode::Integer=0)
            ptr = Ref{Ptr{Cvoid}}()
            status = Ref{Cint}(0)
            ccall(($b,libcfitsio), Cint,
                  (Ref{Ptr{Cvoid}},Ptr{UInt8},Cint,Ref{Cint}),
                  ptr, filename, mode, status)
            fits_assert_ok(status[], filename)
            FITSFile(ptr[])
        end
    end
end

# filename is ignored by the C library
function fits_open_memfile(data::Vector{UInt8}, mode::Integer=0, filename="")
    # Only reading is supported right now
    @assert mode == 0
    ptr = Ref{Ptr{Cvoid}}(C_NULL)
    status = Ref{Cint}(0)
    handle = FITSMemoryHandle(pointer(data),length(data))
    dataptr = Ptr{Ptr{Cvoid}}(pointer_from_objref(handle))
    sizeptr = Ptr{Csize_t}(dataptr+sizeof(Ptr{Cvoid}))
    ccall(("ffomem",libcfitsio), Cint,
      (Ptr{Ptr{Cvoid}},Ptr{UInt8},Cint,Ptr{Ptr{UInt8}},
       Ptr{Csize_t}, Csize_t, Ptr{Cvoid}, Ptr{Cint}),
       ptr, filename, mode, dataptr, sizeptr, 2880, C_NULL, status)
    fits_assert_ok(status[])
    FITSFile(ptr[]), handle
end

"""
    fits_close_file(f::FITSFile)

Close a previously opened FITS file.

"""
function fits_close_file end

"""
    fits_delete_file(f::FITSFile)

Close an opened FITS file (like [`fits_close_file`](@ref)) and removes it from the disk.
"""
function fits_delete_file end

for (a,b) in ((:fits_close_file, "ffclos"),
              (:fits_delete_file,"ffdelt"))
    @eval begin
        function ($a)(f::FITSFile)

            # fits_close_file() is called during garbage collection, but file
            # may already be closed by user, so we need to check if it is open.
            if f.ptr != C_NULL
                status = Ref{Cint}(0)
                ccall(($b,libcfitsio), Cint,
                      (Ptr{Cvoid},Ref{Cint}),
                      f.ptr, status)
                fits_assert_ok(status[])
                f.ptr = C_NULL
            end
        end
    end
end

Base.close(f::FITSFile) = fits_close_file(f)

"""
    fits_file_name(f::FITSFile)

Return the name of the file associated with object `f`.
"""
function fits_file_name(f::FITSFile)
    value = Vector{UInt8}(undef, 1025)
    status = Ref{Cint}(0)
    ccall((:ffflnm,libcfitsio), Cint,
          (Ptr{Cvoid},Ptr{UInt8},Ref{Cint}),
          f.ptr, value, status)
    fits_assert_ok(status[])
    unsafe_string(pointer(value))
end

function fits_file_mode(f::FITSFile)
    result = Ref{Cint}(0)
    status = Ref{Cint}(0)
    ccall(("ffflmd", libcfitsio), Cint, (Ptr{Cvoid}, Ref{Cint}, Ref{Cint}),
          f.ptr, result, status)
    fits_assert_ok(status[])
    result[]
end


# -----------------------------------------------------------------------------
# header access functions

"""
    fits_get_hdrspace(f::FITSFile) -> (keysexist, morekeys)

Return the number of existing keywords (not counting the END keyword)
and the amount of space currently available for more keywords.
"""
function fits_get_hdrspace(f::FITSFile)
    keysexist = Ref{Cint}(0)
    morekeys = Ref{Cint}(0)
    status = Ref{Cint}(0)
    ccall((:ffghsp,libcfitsio), Cint,
        (Ptr{Cvoid},Ref{Cint},Ref{Cint},Ref{Cint}),
        f.ptr, keysexist, morekeys, status)
    fits_assert_ok(status[])
    (keysexist[], morekeys[])
end

function fits_read_key_str(f::FITSFile, keyname::String)
    value = Vector{UInt8}(undef, 71)
    comment = Vector{UInt8}(undef, 71)
    status = Ref{Cint}(0)
    ccall((:ffgkys, libcfitsio), Cint,
          (Ptr{Cvoid}, Ptr{UInt8}, Ptr{UInt8}, Ptr{UInt8}, Ref{Cint}),
          f.ptr, keyname, value, comment, status)
    fits_assert_ok(status[])
    unsafe_string(pointer(value)), unsafe_string(pointer(comment))
end

function fits_read_key_lng(f::FITSFile, keyname::String)
    value = Ref{Clong}(0)
    comment = Vector{UInt8}(undef, 71)
    status = Ref{Cint}(0)
    ccall((:ffgkyj, libcfitsio), Cint,
          (Ptr{Cvoid}, Ptr{UInt8}, Ref{Clong}, Ptr{UInt8}, Ref{Cint}),
          f.ptr, keyname, value, comment, status)
    fits_assert_ok(status[])
    value[], unsafe_string(pointer(comment))
end

function fits_read_keys_lng(f::FITSFile, keyname::String,
                            nstart::Integer, nmax::Integer)
    value = Vector{Clong}(nmax - nstart + 1)
    nfound = Ref{Cint}(0)
    status = Ref{Cint}(0)
    ccall((:ffgknj, libcfitsio), Cint,
          (Ptr{Cvoid}, Ptr{UInt8}, Cint, Cint, Ptr{Clong}, Ref{Cint}, Ref{Cint}),
          f.ptr, keyname, nstart, nmax, value, nfound, status)
    fits_assert_ok(status[])
    value, nfound[]
end

"""
    fits_read_keyword(f::FITSFile, keyname::String) -> (value, comment)

Return the specified keyword.
"""
function fits_read_keyword(f::FITSFile, keyname::String)
    value = Vector{UInt8}(undef, 71)
    comment = Vector{UInt8}(undef, 71)
    status = Ref{Cint}(0)
    ccall((:ffgkey,libcfitsio), Cint,
        (Ptr{Cvoid},Ptr{UInt8},Ptr{UInt8},Ptr{UInt8},Ref{Cint}),
        f.ptr, keyname, value, comment, status)
    fits_assert_ok(status[])
    unsafe_string(pointer(value)), unsafe_string(pointer(comment))
end


"""
    fits_read_record(f::FITSFile, keynum::Int) -> String

Return the nth header record in the CHU. The first keyword in the
header is at `keynum = 1`.
"""
function fits_read_record(f::FITSFile, keynum::Integer)
    card = Vector{UInt8}(undef, 81)
    status = Ref{Cint}(0)
    ccall((:ffgrec,libcfitsio), Cint,
        (Ptr{Cvoid},Cint,Ptr{UInt8},Ref{Cint}),
        f.ptr, keynum, card, status)
    fits_assert_ok(status[])
    unsafe_string(pointer(card))
end


"""
    fits_read_keyn(f::FITSFile, keynum::Int) -> (name, value, comment)

Return the nth header record in the CHU. The first keyword in the header is at `keynum = 1`.
"""
function fits_read_keyn(f::FITSFile, keynum::Integer)
    keyname = Vector{UInt8}(undef, 9)
    value = Vector{UInt8}(undef, 71)
    comment = Vector{UInt8}(undef, 71)
    status = Ref{Cint}(0)
    ccall((:ffgkyn,libcfitsio), Cint,
        (Ptr{Cvoid},Cint,Ptr{UInt8},Ptr{UInt8},Ptr{UInt8},Ref{Cint}),
        f.ptr, keynum, keyname, value, comment, status)
    fits_assert_ok(status[])
    (unsafe_string(pointer(keyname)), unsafe_string(pointer(value)),
     unsafe_string(pointer(comment)))
end

"""
    fits_write_key(f::FITSFile, keyname::String, value, comment::String)

Write a keyword of the appropriate data type into the CHU.
"""
function fits_write_key(f::FITSFile, keyname::String,
                        value::Union{Real,String}, comment::String)
    fits_assert_isascii(keyname)
    fits_assert_isascii(comment)
    cvalue = isa(value,String) ?  value :
             isa(value,Bool) ? Cint[value] : reinterpret(UInt8, [value])
    status = Ref{Cint}(0)
    ccall((:ffpky,libcfitsio), Cint,
        (Ptr{Cvoid},Cint,Ptr{UInt8},Ptr{UInt8},Ptr{UInt8},Ref{Cint}),
        f.ptr, cfitsio_typecode(typeof(value)), keyname,
        cvalue, comment, status)
    fits_assert_ok(status[])
end

function fits_write_date(f::FITSFile)
    status = Ref{Cint}(0)
    ccall((:ffpdat, libcfitsio), Cint, (Ptr{Cvoid}, Ref{Cint}), f.ptr, status)
    fits_assert_ok(status[])
end

function fits_write_comment(f::FITSFile, comment::String)
    fits_assert_isascii(comment)
    status = Ref{Cint}(0)
    ccall((:ffpcom, libcfitsio), Cint, (Ptr{Cvoid}, Ptr{UInt8}, Ref{Cint}),
          f.ptr, comment, status)
    fits_assert_ok(status[])
end

function fits_write_history(f::FITSFile, history::String)
    fits_assert_isascii(history)
    status = Ref{Cint}(0)
    ccall((:ffphis, libcfitsio), Cint, (Ptr{Cvoid}, Ptr{UInt8}, Ref{Cint}),
          f.ptr, history, status)
    fits_assert_ok(status[])
end

# update key: if already present, update it, otherwise add it.
for (a,T,S) in (("ffukys", :String, :(Ptr{UInt8})),
                ("ffukyl", :Bool,        :Cint),
                ("ffukyj", :Integer,     :Int64))
    @eval begin
        function fits_update_key(f::FITSFile, key::String, value::$T,
                                 comment::Union{String, Ptr{Cvoid}}=C_NULL)
            isa(value, String) && fits_assert_isascii(value)
            isa(comment, String) && fits_assert_isascii(comment)
            status = Ref{Cint}(0)
            ccall(($a, libcfitsio), Cint,
                  (Ptr{Cvoid}, Ptr{UInt8}, $S, Ptr{UInt8}, Ref{Cint}),
                  f.ptr, key, value, comment, status)
            fits_assert_ok(status[])
        end
    end
end

function fits_update_key(f::FITSFile, key::String, value::AbstractFloat,
                         comment::Union{String, Ptr{Cvoid}}=C_NULL)
    isa(comment, String) && fits_assert_isascii(comment)
    status = Ref{Cint}(0)
    ccall(("ffukyd", libcfitsio), Cint,
          (Ptr{Cvoid}, Ptr{UInt8}, Cdouble, Cint, Ptr{UInt8}, Ref{Cint}),
          f.ptr, key, value, -15, comment, status)
    fits_assert_ok(status[])
end

function fits_update_key(f::FITSFile, key::String, value::Nothing,
                         comment::Union{String, Ptr{Cvoid}}=C_NULL)
    isa(comment, String) && fits_assert_isascii(comment)
    status = Ref{Cint}(0)
    ccall(("ffukyu", libcfitsio), Cint,
          (Ptr{Cvoid}, Ptr{UInt8}, Ptr{UInt8}, Ref{Cint}),
          f.ptr, key, comment, status)
    fits_assert_ok(status[])
end

"""
    fits_write_record(f::FITSFile, card::String)

Write a user specified keyword record into the CHU.
"""
function fits_write_record(f::FITSFile, card::String)
    fits_assert_isascii(card)
    status = Ref{Cint}(0)
    ccall((:ffprec,libcfitsio), Cint,
        (Ptr{Cvoid},Ptr{UInt8},Ref{Cint}),
        f.ptr, card, status)
    fits_assert_ok(status[])
end

"""
    fits_delete_record(f::FITSFile, keynum::Int)

Delete the keyword record at the specified index.
"""
function fits_delete_record(f::FITSFile, keynum::Integer)
    status = Ref{Cint}(0)
    ccall((:ffdrec,libcfitsio), Cint,
        (Ptr{Cvoid},Cint,Ref{Cint}),
        f.ptr, keynum, status)
    fits_assert_ok(status[])
end

"""
    fits_delete_key(f::FITSFile, keyname::String)

Delete the keyword named `keyname`.
"""
function fits_delete_key(f::FITSFile, keyname::String)
    status = Ref{Cint}(0)
    ccall((:ffdkey,libcfitsio), Cint,
        (Ptr{Cvoid},Ptr{UInt8},Ref{Cint}),
        f.ptr, keyname, status)
    fits_assert_ok(status[])
end

"""
    fits_hdr2str(f::FITSFile, nocomments::Bool=false)

Return the header of the CHDU as a string. If `nocomments` is `true`, comment
cards are stripped from the output.
"""
function fits_hdr2str(f::FITSFile, nocomments::Bool=false)
    status = Ref{Cint}(0)
    header = Ref{Ptr{UInt8}}()
    nkeys = Ref{Cint}(0)
    ccall((:ffhdr2str, libcfitsio), Cint,
          (Ptr{Cvoid}, Cint, Ref{Ptr{UInt8}}, Cint,
           Ptr{Ptr{UInt8}}, Ref{Cint}, Ref{Cint}),
          f.ptr, nocomments, C_NULL, 0, header, nkeys, status)
    result = unsafe_string(header[])

    # free header pointer allocated by cfitsio (result is a copy)
    ccall((:fffree, libcfitsio), Ref{Cint},
          (Ptr{UInt8}, Ref{Cint}),
          header[], status)
    fits_assert_ok(status[])
    result
end


# -----------------------------------------------------------------------------
# HDU info functions and moving the current HDU

function hdu_int_to_type(hdu_type_int)
    if hdu_type_int == 0
        return :image_hdu
    elseif hdu_type_int == 1
        return :ascii_table
    elseif hdu_type_int == 2
        return :binary_table
    end

    :unknown
end

"""
    fits_movabs_hdu(f::FITSFile, hduNum::Integer)

Change the current HDU to the value specified by `hduNum`, and return a symbol
describing the type of the HDU.

Possible symbols are: `image_hdu`, `ascii_table`, or `binary_table`.
The value of `hduNum` must range between 1 and the value returned by
[`fits_get_num_hdus`](@ref).
"""
function fits_movabs_hdu end

"""
    fits_movrel_hdu(f::FITSFile, hduNum::Integer)

Change the current HDU by moving forward or backward by `hduNum` HDUs
(positive means forward), and return the same as [`fits_movabs_hdu`](@ref).
"""
function fits_movrel_hdu end
for (a,b) in ((:fits_movabs_hdu,"ffmahd"),
              (:fits_movrel_hdu,"ffmrhd"))
    @eval begin
        function ($a)(f::FITSFile, hduNum::Integer)
            hdu_type = Ref{Cint}(0)
            status = Ref{Cint}(0)
            ccall(($b,libcfitsio), Cint,
                  (Ptr{Cvoid}, Cint, Ref{Cint}, Ref{Cint}),
                  f.ptr, hduNum, hdu_type, status)
            fits_assert_ok(status[])
            hdu_int_to_type(hdu_type[])
        end
    end
end

"""
    fits_movnam_hdu(f::FITSFile, extname::String, extver::Integer=0,
                    hdu_type_int::Integer=-1)

Change the current HDU by moving to the (first) HDU which has the specified
extension type and EXTNAME and EXTVER keyword values (or HDUNAME and HDUVER keywords).

If `extver` is 0 (the default) then the EXTVER keyword is ignored and the first HDU
with a matching EXTNAME (or HDUNAME) keyword will be found. If `hdu_type_int`
is -1 (the default) only the extname and extver values will be used to locate the
correct extension. If no matching HDU is found in the file, the current HDU will
remain unchanged.
"""
function fits_movnam_hdu(f::FITSFile, extname::String, extver::Integer=0,
                         hdu_type::Integer=-1)
    status = Ref{Cint}(0)
    ccall((:ffmnhd,libcfitsio), Cint,
          (Ptr{Cvoid}, Cint, Ptr{UInt8}, Cint, Ref{Cint}),
          f.ptr, hdu_type, extname, extver, status)
    fits_assert_ok(status[])
end

function fits_get_hdu_num(f::FITSFile)
    hdunum = Ref{Cint}(0)
    ccall((:ffghdn,libcfitsio), Cint,
          (Ptr{Cvoid}, Ref{Cint}),
          f.ptr, hdunum)
    hdunum[]
end

function fits_get_hdu_type(f::FITSFile)
    hdutype = Ref{Cint}(0)
    status = Ref{Cint}(0)
    ccall((:ffghdt, libcfitsio), Cint,
          (Ptr{Cvoid}, Ref{Cint}, Ref{Cint}),
          f.ptr, hdutype, status)
    fits_assert_ok(status[])
    hdu_int_to_type(hdutype[])
end

# -----------------------------------------------------------------------------
# image HDU functions

"""
    fits_get_img_size(f::FITSFile)

Get the dimensions of the image.
"""
function fits_get_img_size end
for (a, b) in ((:fits_get_img_type,      "ffgidt"),
               (:fits_get_img_equivtype, "ffgiet"),
               (:fits_get_img_dim,       "ffgidm"))
    @eval function ($a)(f::FITSFile)
        result = Ref{Cint}(0)
        status = Ref{Cint}(0)
        ccall(($b, libcfitsio), Cint, (Ptr{Cvoid}, Ref{Cint}, Ref{Cint}),
              f.ptr, result, status)
        fits_assert_ok(status[])
        result[]
    end
end

"""
    fits_create_img(f::FITSFile, t::Type, naxes::Vector{Int})

Create a new primary array or IMAGE extension with a specified data type and size.
"""
function fits_create_img(f::FITSFile, ::Type{T},
                         naxes::Vector{S}) where {T, S<:Integer}
    status = Ref{Cint}(0)
    ccall((:ffcrimll, libcfitsio), Cint,
          (Ptr{Cvoid}, Cint, Cint, Ptr{Int64}, Ref{Cint}),
          f.ptr, bitpix_from_type(T), length(naxes),
          Vector{Int64}(naxes), status)
    fits_assert_ok(status[])
end

"""
    fits_write_pix(f::FITSFile, fpixel::Vector{Int}, nelements::Int, data::Array)

Write pixels from `data` into the FITS file.
"""
function fits_write_pix(f::FITSFile, fpixel::Vector{S},
                        nelements::Integer, data::Array{T}) where {S<:Integer,T}
    status = Ref{Cint}(0)
    ccall((:ffppxll, libcfitsio), Cint,
          (Ptr{Cvoid}, Cint, Ptr{Int64}, Int64, Ptr{Cvoid}, Ref{Cint}),
          f.ptr, cfitsio_typecode(T), Vector{Int64}(fpixel),
          nelements, data, status)
    fits_assert_ok(status[])
end

function fits_write_pix(f::FITSFile, data::Array)
    fits_write_pix(f, ones(Int64, length(size(data))), length(data), data)
end

function fits_read_pix(f::FITSFile, fpixel::Vector{S},
                       nelements::Int, nullval::T,
                       data::Array{T}) where {S<:Integer,T}
    anynull = Ref{Cint}(0)
    status = Ref{Cint}(0)
    ccall((:ffgpxvll, libcfitsio), Cint,
          (Ptr{Cvoid}, Cint, Ptr{Int64}, Int64, Ref{Cvoid}, Ptr{Cvoid},
           Ref{Cint}, Ref{Cint}),
          f.ptr, cfitsio_typecode(T), Vector{Int64}(fpixel),
          nelements, nullval, data, anynull, status)
    fits_assert_ok(status[])
    anynull[]
end

"""
    fits_read_pix(f::FITSFile, fpixel::Vector{Int}, nelements::Int, data::Array)

Read pixels from the FITS file into `data`.
"""
function fits_read_pix(f::FITSFile, fpixel::Vector{S},
                       nelements::Int, data::Array{T}) where {S<:Integer,T}
    anynull = Ref{Cint}(0)
    status = Ref{Cint}(0)
    ccall((:ffgpxvll, libcfitsio), Cint,
          (Ptr{Cvoid}, Cint, Ptr{Int64}, Int64, Ptr{Cvoid}, Ptr{Cvoid},
           Ref{Cint}, Ref{Cint}),
          f.ptr, cfitsio_typecode(T), Vector{Int64}(fpixel),
          nelements, C_NULL, data, anynull, status)
    fits_assert_ok(status[])
    anynull[]
end

function fits_read_pix(f::FITSFile, data::Array)
    fits_read_pix(f, ones(Int64,length(size(data))), length(data), data)
end

function fits_read_subset(
             f::FITSFile, fpixel::Vector{S1}, lpixel::Vector{S2},
             inc::Vector{S3}, data::Array{T}) where {S1<:Integer,S2<:Integer,S3<:Integer,T}
    anynull = Ref{Cint}(0)
    status = Ref{Cint}(0)
    ccall((:ffgsv, libcfitsio), Cint,
          (Ptr{Cvoid}, Cint, Ptr{Clong}, Ptr{Clong}, Ptr{Clong}, Ptr{Cvoid},
           Ptr{Cvoid}, Ref{Cint}, Ref{Cint}),
          f.ptr, cfitsio_typecode(T),
          Vector{Clong}(fpixel),
          Vector{Clong}(lpixel),
          Vector{Clong}(inc),
          C_NULL, data, anynull, status)
    fits_assert_ok(status[])
    anynull[]
end

function fits_copy_image_section(fin::FITSFile, fout::FITSFile,
                                 section::String)
    status = Ref{Cint}(0)
    ccall((:fits_copy_image_section,libcfitsio), Cint,
          (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{UInt8}, Ref{Cint}),
          fin.ptr, fout.ptr, section, status)
    fits_assert_ok(status[])
end


# -----------------------------------------------------------------------------
# ASCII/binary table HDU functions

# The three fields are: ttype, tform, tunit (CFITSIO's terminology)
const ColumnDef = Tuple{String, String, String}

"""
    fits_create_binary_tbl(f::FITSFile, numrows::Integer, coldefs::Array{ColumnDef},
                           extname::String)

Append a new HDU containing a binary table. The meaning of the parameters is the same
as in a call to [`fits_create_ascii_tbl`](@ref).

In general, one should pick this function for creating tables in a new HDU,
as binary tables require less space on the disk and are more efficient to read and write.
(Moreover, a few datatypes are not supported in ASCII tables).
"""
function fits_create_binary_tbl end

"""
    fits_create_ascii_tbl(f::FITSFile, numrows::Integer, coldefs::Array{ColumnDef},
                          extname::String)

Append a new HDU containing an ASCII table.

The table will have `numrows` rows (this parameter can be set to zero), each
initialized with the default value. In order to create a table, the programmer
must specify the characteristics of each column. The columns are specified by the
`coldefs` variable, which is an array of tuples.
Each tuple must have three string fields:

1. The name of the column.
2. The data type and the repetition count. It must be a string made by a number
   (the repetition count) followed by a letter specifying the type (in the example
   above, `D` stands for `Float64`, `E` stands for `Float32`, `A` stands for `Char`).
   Refer to the CFITSIO documentation for more information about the syntax of this
   parameter.
3. The measure unit of this field. This is used only as a comment.

The value of `extname` sets the "extended name" of the table, i.e., a string
that in some situations can be used to refer to the HDU itself.

Note that, unlike for binary tables, CFITSIO puts some limitations to the
types that can be used in an ASCII table column. Refer to the CFITSIO manual
for further information.

See also [`fits_create_binary_tbl`](@ref) for a similar function which
creates binary tables.
"""
function fits_create_ascii_tbl end
for (a,b) in ((:fits_create_binary_tbl, 2),
              (:fits_create_ascii_tbl,  1))
    @eval begin
        function ($a)(f::FITSFile, numrows::Integer,
                      coldefs::Array{ColumnDef}, extname::String)

            # Ensure that extension name, column names and units are
            # ASCII, as these get written to the file. We don't check
            # need to check that tform is ASCII because presumably
            # cfitsio will thrown an appropriate error if it doesn't
            # recognize the tform string.
            fits_assert_isascii(extname)
            for coldef in coldefs
                fits_assert_isascii(coldef[1])
                fits_assert_isascii(coldef[3])
            end

            # get length and convert coldefs to three arrays of Ptr{Uint8}
            ntype = length(coldefs)
            ttype = [pointer(x[1]) for x in coldefs]
            tform = [pointer(x[2]) for x in coldefs]
            tunit = [pointer(x[3]) for x in coldefs]
            status = Ref{Cint}(0)

            ccall(("ffcrtb", libcfitsio), Cint,
                  (Ptr{Cvoid}, Cint, Int64, Cint,
                   Ptr{Ptr{UInt8}}, Ptr{Ptr{UInt8}},
                   Ptr{Ptr{UInt8}}, Ptr{UInt8}, Ref{Cint}),
                  f.ptr, $b, numrows, ntype,
                  ttype, tform, tunit, extname,
                  status)
            fits_assert_ok(status[])
        end
    end
end

"""
    fits_get_num_hdus(f::FITSFile)

Return the number of HDUs in the file.
"""
function fits_get_num_hdus end
for (a,b,T) in ((:fits_get_num_cols,  "ffgncl",  :Cint),
                (:fits_get_num_hdus,  "ffthdu",  :Cint),
                (:fits_get_rowsize,   "ffgrsz",  :Clong))
    @eval begin
        function ($a)(f::FITSFile)
            result = Ref{$T}(0)
            status = Ref{Cint}(0)
            ccall(($b,libcfitsio), Cint,
                  (Ptr{Cvoid}, Ref{$T}, Ref{Cint}),
                  f.ptr, result, status)
            fits_assert_ok(status[])
            result[]
        end
    end
end

function fits_get_colnum(f::FITSFile, tmplt::String; case_sensitive::Bool=true)
    result = Ref{Cint}(0)
    status = Ref{Cint}(0)

    # Second argument is case-sensitivity of search: 0 = case-insensitive
    #                                                1 = case-sensitive
    ccall(("ffgcno", libcfitsio), Cint,
          (Ptr{Cvoid}, Cint, Ptr{UInt8}, Ref{Cint}, Ref{Cint}),
          f.ptr, case_sensitive, tmplt, result, status)
    fits_assert_ok(status[])
    return result[]
end

# The following block are all functions that have separate variants for Clong
# and 64-bit integers in cfitsio. Rather than providing both of these, we
# provide only one according to the native integer type on the platform.
if  promote_type(Int, Clong) == Clong
    T = Clong
    ffgtdm = "ffgtdm"
    ffgnrw = "ffgnrw"
    ffptdm = "ffptdm"
    ffgtcl = "ffgtcl"
    ffeqty = "ffeqty"
    ffgdes = "ffgdes"
    ffgisz = "ffgisz"
else
    T = Int64
    ffgtdm = "ffgtdmll"
    ffgnrw = "ffgnrwll"
    ffptdm = "ffptdmll"
    ffgtcl = "ffgtclll"
    ffeqty = "ffeqtyll"
    ffgdes = "ffgdesll"
    ffgisz = "ffgiszll"
end

"""
    fits_get_coltype(f::FITSFile, colnum::Integer)

Provided that the current HDU contains either an ASCII or binary table, return
information about the column at position `colnum` (counting from 1).

Return is a tuple containing

- `typecode`: CFITSIO integer type code of the column.
- `repcount`: Repetition count for the column.
- `width`: Width of an individual element.
"""
function fits_get_coltype end

@eval begin
    function fits_get_coltype(ff::FITSFile, colnum::Integer)
        typecode = Ref{Cint}(0)
        repcnt = Ref{$T}(0)
        width = Ref{$T}(0)
        status = Ref{Cint}(0)
        ccall(($ffgtcl,libcfitsio), Cint,
              (Ptr{Cvoid}, Cint, Ref{Cint}, Ref{$T}, Ref{$T}, Ref{Cint}),
              ff.ptr, colnum, typecode, repcnt, width, status)
        fits_assert_ok(status[])
        return Int(typecode[]), Int(repcnt[]), Int(width[])
    end

    function fits_get_eqcoltype(ff::FITSFile, colnum::Integer)
        typecode = Ref{Cint}(0)
        repcnt = Ref{$T}(0)
        width = Ref{$T}(0)
        status = Ref{Cint}(0)
        ccall(($ffeqty,libcfitsio), Cint,
              (Ptr{Cvoid}, Cint, Ref{Cint}, Ref{$T}, Ref{$T}, Ref{Cint}),
              ff.ptr, colnum, typecode, repcnt, width, status)
        fits_assert_ok(status[])
        return Int(typecode[]), Int(repcnt[]), Int(width[])
    end

    function fits_get_img_size(f::FITSFile)
        ndim = fits_get_img_dim(f)
        naxes = Vector{$T}(ndim)
        status = Ref{Cint}(0)
        ccall(($ffgisz, libcfitsio), Cint,
              (Ptr{Cvoid}, Cint, Ptr{$T}, Ref{Cint}),
              f.ptr, ndim, naxes, status)
        fits_assert_ok(status[])
        naxes
    end

    function fits_get_num_rows(f::FITSFile)
        result = Ref{$T}(0)
        status = Ref{Cint}(0)
        ccall(($ffgnrw, libcfitsio), Cint,
              (Ptr{Cvoid}, Ref{$T}, Ref{Cint}),
              f.ptr, result, status)
        fits_assert_ok(status[])
        return Int(result[])
    end

    # `fits_read_tdim` returns the dimensions of a table column in a
    # binary table. Normally this information is given by the TDIMn
    # keyword, but if this keyword is not present then this routine
    # returns `[r]` with `r` equals to the repeat count in the TFORM
    # keyword.
    function fits_read_tdim(ff::FITSFile, colnum::Integer)
        naxes = Vector{$T}(undef, 99)  # 99 is the maximum allowed number of axes
        naxis = Ref{Cint}(0)
        status = Ref{Cint}(0)
        ccall(($ffgtdm,libcfitsio), Cint,
              (Ptr{Cvoid}, Cint, Cint, Ref{Cint}, Ptr{$T}, Ref{Cint}),
              ff.ptr, colnum, length(naxes), naxis, naxes, status)
        fits_assert_ok(status[])
        return naxes[1:naxis[]]
    end

    function fits_write_tdim(ff::FITSFile, colnum::Integer,
                                 naxes::Array{$T})
        status = Ref{Cint}(0)
        ccall(($ffptdm, libcfitsio), Cint,
              (Ptr{Cvoid}, Cint, Cint, Ptr{$T}, Ref{Cint}),
              ff.ptr, colnum, length(naxes), naxes, status)
        fits_assert_ok(status[])
    end

    function fits_read_descript(f::FITSFile, colnum::Integer, rownum::Integer)
        repeat = Ref{$T}(0)
        offset = Ref{$T}(0)
        status = Ref{Cint}(0)
        ccall(($ffgdes, libcfitsio), Cint,
              (Ptr{Cvoid}, Cint, Int64, Ref{$T}, Ref{$T}, Ref{Cint}),
              f.ptr, colnum, rownum, repeat, offset, status)
        fits_assert_ok(status[])
        return Int(repeat[]), Int(offset[])
    end
end

"""
    fits_read_col(f, colnum, firstrow, firstelem, data)

Read data from one column of an ASCII/binary table and convert the data into the
specified type `T`.

### Arguments ###

* `f::FITSFile`: the file to be read.
* `colnum::Integer`: the column number, where the value of the first column is `1`.
* `firstrow::Integer`: the elements to be read start from this row.
* `firstelem::Integer`: specifies which is the first element to be read, when each
  cell contains more than one element (i.e., the "repetition count" of the field is
  greater than one).
* `data::Array`: at the end of the call, this will be filled with the elements read
  from the column. The length of the array gives the overall number of elements.
"""
function fits_read_col(f::FITSFile,
                       colnum::Integer,
                       firstrow::Integer,
                       firstelem::Integer,
                       data::Array{String})

    # get width: number of characters in each string
    typecode, repcount, width = fits_get_eqcoltype(f, colnum)

    # ensure that data are strings, otherwise cfitsio will try to write
    # formatted strings, which have widths given by fits_get_col_display_width
    # not by the repeat value from fits_get_coltype.
    abs(typecode) == 16 || error("not a string column")

    # create an array of character buffers of the correct width
    buffers = [Vector{UInt8}(width) for i in 1:length(data)]

    # Call the CFITSIO function
    anynull = Ref{Cint}(0)
    status = Ref{Cint}(0)
    ccall((:ffgcvs, libcfitsio), Cint,
          (Ptr{Cvoid}, Cint, Int64, Int64, Int64,
           Ptr{UInt8}, Ptr{Ptr{UInt8}}, Ref{Cint}, Ref{Cint}),
          f.ptr, colnum, firstrow, firstelem, length(data),
          " ", buffers, anynull, status)
    fits_assert_ok(status[])

    # Create strings out of the buffers, terminating at null characters.
    # Note that `String(x)` does not copy the buffer x.
    for i in 1:length(data)
        zeropos = search(buffers[i], 0x00)
        data[i] = (zeropos >= 1) ? String(buffers[i][1:(zeropos-1)]) :
                                   String(buffers[i])
    end
end

function fits_read_col(f::FITSFile,
                       colnum::Integer,
                       firstrow::Integer,
                       firstelem::Integer,
                       data::Array{T}) where T
    anynull = Ref{Cint}(0)
    status = Ref{Cint}(0)
    ccall((:ffgcv,libcfitsio), Cint,
          (Ptr{Cvoid}, Cint, Cint, Int64, Int64, Int64,
           Ptr{Cvoid}, Ptr{Cvoid}, Ref{Cint}, Ref{Cint}),
          f.ptr, cfitsio_typecode(T), colnum,
          firstrow, firstelem, length(data),
          C_NULL, data, anynull, status)
    fits_assert_ok(status[])
end

"""
    fits_write_col(f, colnum, firstrow, firstelem, data)

Write some data in one column of a ASCII/binary table.

If there is no room for the elements, new rows will be created. (It is therefore
useless to call [`fits_insert_rows`](@ref) if you only need to *append* elements
to the end of a table.)

* `f::FITSFile`: the file in which data will be written.
* `colnum::Integer`: the column number, where the value of the first column is `1`.
* `firstrow::Integer`: the data wil be written from this row onwards.
* `firstelem::Integer`: specifies the position in the row where the first element
  will be written.
* `data::Array`: contains the elements that are to be written to the column of the table.
"""
function fits_write_col(f::FITSFile,
                        colnum::Integer,
                        firstrow::Integer,
                        firstelem::Integer,
                        data::Array{String})
    for el in data; fits_assert_isascii(el); end
    status = Ref{Cint}(0)
    ccall((:ffpcls, libcfitsio), Cint,
          (Ptr{Cvoid}, Cint, Int64, Int64, Int64,
           Ptr{Ptr{UInt8}}, Ref{Cint}),
          f.ptr, colnum, firstrow, firstelem, length(data),
          data, status)
    fits_assert_ok(status[])
end

function fits_write_col(f::FITSFile,
                        colnum::Integer,
                        firstrow::Integer,
                        firstelem::Integer,
                        data::Array{T}) where T
    status = Ref{Cint}(0)
    ccall((:ffpcl, libcfitsio), Cint,
          (Ptr{Cvoid}, Cint, Cint, Int64, Int64, Int64,
           Ptr{Cvoid}, Ref{Cint}),
          f.ptr, cfitsio_typecode(T), colnum,
          firstrow, firstelem, length(data),
          data, status)
    fits_assert_ok(status[])
end

"""
    fits_insert_rows(f::FITSFile, firstrow::Integer, nrows::Integer)

Insert a number of rows equal to `nrows` after the row number `firstrow`.

The elements in each row are initialized to their default value: you can
modify them later using [`fits_write_col`](@ref).

Since the first row is at position 1, in order to insert rows *before*
the first one `firstrow` must be equal to zero.
"""
function fits_insert_rows end

"""
    fits_delete_rows(f::FITSFile, firstrow::integer, nrows::Integer)

Delete `nrows` rows, starting from the one at position `firstrow`. The index of
the first row is 1.
"""
function fits_delete_rows end

for (a,b) in ((:fits_insert_rows, "ffirow"),
              (:fits_delete_rows, "ffdrow"))
    @eval begin
        function ($a)(f::FITSFile, firstrow::Integer, nrows::Integer)
            status = Ref{Cint}(0)
            ccall(($b,libcfitsio), Cint,
                  (Ptr{Cvoid}, Int64, Int64, Ref{Cint}),
                  f.ptr, firstrow, nrows, status)
            fits_assert_ok(status[])
        end
    end
end

end # module
