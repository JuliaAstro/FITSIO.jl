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
#       16  TSTRING         Compat.ASCIIString
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

isdefined(Base, :__precompile__) && __precompile__()

module Libcfitsio

using Compat

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
                  (Compat.ASCIIString, 16),
                  (Cushort,     20),
                  (Cshort,      21),
                  (Cuint,       30),
                  (Cint,        31),
                  (Int64,       81),
                  (Float32,     42),
                  (Float64,     82),
                  (Complex64,   83),
                  (Complex128, 163))
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

type FITSFile
    ptr::Ptr{Void}

    function FITSFile(ptr::Ptr{Void})
        f = new(ptr)
        finalizer(f, fits_close_file)
        f
    end
end

# FITS wants to be able to update the ptr, so keep them
# in a mutable struct
type FITSMemoryHandle
    ptr::Ptr{Void}
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
    msg = Vector{UInt8}(31)
    ccall((:ffgerr,libcfitsio), Void, (Cint,Ptr{UInt8}), status, msg)
    unsafe_string(pointer(msg))
end

function fits_assert_ok(status::Cint)
    if status != 0
        error(fits_get_errstatus(status))
    end
end

fits_get_version() = ccall((:ffvers, libcfitsio), Cfloat, (Ptr{Cfloat},), &0.)

# -----------------------------------------------------------------------------
# file access & info functions

function fits_create_file(filename::AbstractString)
    ptr = Ref{Ptr{Void}}()
    status = Ref{Cint}(0)
    ccall((:ffinit,libcfitsio), Cint, (Ref{Ptr{Void}},Ptr{UInt8},Ref{Cint}),
          ptr, filename, status)
    fits_assert_ok(status[])
    FITSFile(ptr[])
end

fits_clobber_file(filename::AbstractString) = fits_create_file("!"*filename)

for (a,b) in ((:fits_open_data, "ffdopn"),
              (:fits_open_file, "ffopen"),
              (:fits_open_image,"ffiopn"),
              (:fits_open_table,"fftopn"))
    @eval begin
        function ($a)(filename::AbstractString, mode::Integer=0)
            ptr = Ref{Ptr{Void}}()
            status = Ref{Cint}(0)
            ccall(($b,libcfitsio), Cint,
                  (Ref{Ptr{Void}},Ptr{UInt8},Cint,Ref{Cint}),
                  ptr, filename, mode, status)
            fits_assert_ok(status[])
            FITSFile(ptr[])
        end
    end
end

# filename is ignored by the C library
function fits_open_memfile(data::Vector{UInt8}, mode::Integer=0, filename="")
    # Only reading is supported right now
    @assert mode == 0
    ptr = Ref{Ptr{Void}}(C_NULL)
    status = Ref{Cint}(0)
    handle = FITSMemoryHandle(pointer(data),length(data))
    dataptr = Ptr{Ptr{Void}}(pointer_from_objref(handle))
    sizeptr = Ptr{Csize_t}(dataptr+sizeof(Ptr{Void}))
    ccall(("ffomem",libcfitsio), Cint,
      (Ptr{Ptr{Void}},Ptr{UInt8},Cint,Ptr{Ptr{UInt8}},
       Ptr{Csize_t}, Csize_t, Ptr{Void}, Ptr{Cint}),
       ptr, filename, mode, dataptr, sizeptr, 2880, C_NULL, status)
    fits_assert_ok(status[])
    FITSFile(ptr[]), handle
end

for (a,b) in ((:fits_close_file, "ffclos"),
              (:fits_delete_file,"ffdelt"))
    @eval begin
        function ($a)(f::FITSFile)

            # fits_close_file() is called during garbage collection, but file
            # may already be closed by user, so we need to check if it is open.
            if f.ptr != C_NULL
                status = Ref{Cint}(0)
                ccall(($b,libcfitsio), Cint,
                      (Ptr{Void},Ref{Cint}),
                      f.ptr, status)
                fits_assert_ok(status[])
                f.ptr = C_NULL
            end
        end
    end
end

close(f::FITSFile) = fits_close_file(f)

function fits_file_name(f::FITSFile)
    value = Vector{UInt8}(1025)
    status = Ref{Cint}(0)
    ccall((:ffflnm,libcfitsio), Cint,
          (Ptr{Void},Ptr{UInt8},Ref{Cint}),
          f.ptr, value, status)
    fits_assert_ok(status[])
    unsafe_string(pointer(value))
end

function fits_file_mode(f::FITSFile)
    result = Ref{Cint}(0)
    status = Ref{Cint}(0)
    ccall(("ffflmd", libcfitsio), Cint, (Ptr{Void}, Ref{Cint}, Ref{Cint}),
          f.ptr, result, status)
    fits_assert_ok(status[])
    result[]
end


# -----------------------------------------------------------------------------
# header access functions

function fits_get_hdrspace(f::FITSFile)
    keysexist = Ref{Cint}(0)
    morekeys = Ref{Cint}(0)
    status = Ref{Cint}(0)
    ccall((:ffghsp,libcfitsio), Cint,
        (Ptr{Void},Ref{Cint},Ref{Cint},Ref{Cint}),
        f.ptr, keysexist, morekeys, status)
    fits_assert_ok(status[])
    (keysexist[], morekeys[])
end

function fits_read_key_str(f::FITSFile, keyname::Compat.ASCIIString)
    value = Vector{UInt8}(71)
    comment = Vector{UInt8}(71)
    status = Ref{Cint}(0)
    ccall((:ffgkys, libcfitsio), Cint,
          (Ptr{Void}, Ptr{UInt8}, Ptr{UInt8}, Ptr{UInt8}, Ref{Cint}),
          f.ptr, keyname, value, comment, status)
    fits_assert_ok(status[])
    unsafe_string(pointer(value)), unsafe_string(pointer(comment))
end

function fits_read_key_lng(f::FITSFile, keyname::Compat.ASCIIString)
    value = Ref{Clong}(0)
    comment = Vector{UInt8}(71)
    status = Ref{Cint}(0)
    ccall((:ffgkyj, libcfitsio), Cint,
          (Ptr{Void}, Ptr{UInt8}, Ref{Clong}, Ptr{UInt8}, Ref{Cint}),
          f.ptr, keyname, value, comment, status)
    fits_assert_ok(status[])
    value[], unsafe_string(pointer(comment))
end

function fits_read_keys_lng(f::FITSFile, keyname::Compat.ASCIIString,
                            nstart::Int, nmax::Int)
    value = Vector{Clong}(nmax - nstart + 1)
    nfound = Ref{Cint}(0)
    status = Ref{Cint}(0)
    ccall((:ffgknj, libcfitsio), Cint,
          (Ptr{Void}, Ptr{UInt8}, Cint, Cint, Ptr{Clong}, Ref{Cint}, Ref{Cint}),
          f.ptr, keyname, nstart, nmax, value, nfound, status)
    fits_assert_ok(status[])
    value, nfound[]
end

function fits_read_keyword(f::FITSFile, keyname::Compat.ASCIIString)
    value = Vector{UInt8}(71)
    comment = Vector{UInt8}(71)
    status = Ref{Cint}(0)
    ccall((:ffgkey,libcfitsio), Cint,
        (Ptr{Void},Ptr{UInt8},Ptr{UInt8},Ptr{UInt8},Ref{Cint}),
        f.ptr, keyname, value, comment, status)
    fits_assert_ok(status[])
    unsafe_string(pointer(value)), unsafe_string(pointer(comment))
end

function fits_read_record(f::FITSFile, keynum::Integer)
    card = Vector{UInt8}(81)
    status = Ref{Cint}(0)
    ccall((:ffgrec,libcfitsio), Cint,
        (Ptr{Void},Cint,Ptr{UInt8},Ref{Cint}),
        f.ptr, keynum, card, status)
    fits_assert_ok(status[])
    unsafe_string(pointer(card))
end

function fits_read_keyn(f::FITSFile, keynum::Integer)
    keyname = Vector{UInt8}(9)
    value = Vector{UInt8}(71)
    comment = Vector{UInt8}(71)
    status = Ref{Cint}(0)
    ccall((:ffgkyn,libcfitsio), Cint,
        (Ptr{Void},Cint,Ptr{UInt8},Ptr{UInt8},Ptr{UInt8},Ref{Cint}),
        f.ptr, keynum, keyname, value, comment, status)
    fits_assert_ok(status[])
    (unsafe_string(pointer(keyname)), unsafe_string(pointer(value)),
     unsafe_string(pointer(comment)))
end

function fits_write_key(f::FITSFile, keyname::Compat.ASCIIString,
                        value::Union{Real,Compat.ASCIIString},
                        comment::Compat.ASCIIString)
    cvalue = isa(value,Compat.ASCIIString) ?  value :
             isa(value,Bool) ? Cint[value] : reinterpret(UInt8, [value])
    status = Ref{Cint}(0)
    ccall((:ffpky,libcfitsio), Cint,
        (Ptr{Void},Cint,Ptr{UInt8},Ptr{UInt8},Ptr{UInt8},Ref{Cint}),
        f.ptr, cfitsio_typecode(typeof(value)), keyname,
        cvalue, comment, status)
    fits_assert_ok(status[])
end

function fits_write_date(f::FITSFile)
    status = Ref{Cint}(0)
    ccall((:ffpdat, libcfitsio), Cint, (Ptr{Void}, Ref{Cint}), f.ptr, status)
    fits_assert_ok(status[])
end

function fits_write_comment(f::FITSFile, comment::Compat.ASCIIString)
    status = Ref{Cint}(0)
    ccall((:ffpcom, libcfitsio), Cint, (Ptr{Void}, Ptr{UInt8}, Ref{Cint}),
          f.ptr, comment, status)
    fits_assert_ok(status[])
end

function fits_write_history(f::FITSFile, history::Compat.ASCIIString)
    status = Ref{Cint}(0)
    ccall((:ffphis, libcfitsio), Cint, (Ptr{Void}, Ptr{UInt8}, Ref{Cint}),
          f.ptr, history, status)
    fits_assert_ok(status[])
end

# update key: if already present, update it, otherwise add it.
for (a,T,S) in (("ffukys", :(Compat.ASCIIString), :(Ptr{UInt8})),
                ("ffukyl", :Bool,        :Cint),
                ("ffukyj", :Integer,     :Int64))
    @eval begin
        function fits_update_key(f::FITSFile, key::Compat.ASCIIString, value::$T,
                                 comment::Union{Compat.ASCIIString, Ptr{Void}}=C_NULL)
            status = Ref{Cint}(0)
            ccall(($a, libcfitsio), Cint,
                  (Ptr{Void}, Ptr{UInt8}, $S, Ptr{UInt8}, Ref{Cint}),
                  f.ptr, key, value, comment, status)
            fits_assert_ok(status[])
        end
    end
end

function fits_update_key(f::FITSFile, key::Compat.ASCIIString, value::AbstractFloat,
                         comment::Union{Compat.ASCIIString, Ptr{Void}}=C_NULL)
    status = Ref{Cint}(0)
    ccall(("ffukyd", libcfitsio), Cint,
          (Ptr{Void}, Ptr{UInt8}, Cdouble, Cint, Ptr{UInt8}, Ref{Cint}),
          f.ptr, key, value, -15, comment, status)
    fits_assert_ok(status[])
end

function fits_update_key(f::FITSFile, key::Compat.ASCIIString, value::Void,
                         comment::Union{Compat.ASCIIString, Ptr{Void}}=C_NULL)
    status = Ref{Cint}(0)
    ccall(("ffukyu", libcfitsio), Cint,
          (Ptr{Void}, Ptr{UInt8}, Ptr{UInt8}, Ref{Cint}),
          f.ptr, key, comment, status)
    fits_assert_ok(status[])
end

function fits_write_record(f::FITSFile, card::Compat.ASCIIString)
    status = Ref{Cint}(0)
    ccall((:ffprec,libcfitsio), Cint,
        (Ptr{Void},Ptr{UInt8},Ref{Cint}),
        f.ptr, card, status)
    fits_assert_ok(status[])
end

function fits_delete_record(f::FITSFile, keynum::Integer)
    status = Ref{Cint}(0)
    ccall((:ffdrec,libcfitsio), Cint,
        (Ptr{Void},Cint,Ref{Cint}),
        f.ptr, keynum, status)
    fits_assert_ok(status[])
end

function fits_delete_key(f::FITSFile, keyname::Compat.ASCIIString)
    status = Ref{Cint}(0)
    ccall((:ffdkey,libcfitsio), Cint,
        (Ptr{Void},Ptr{UInt8},Ref{Cint}),
        f.ptr, keyname, status)
    fits_assert_ok(status[])
end

function fits_hdr2str(f::FITSFile, nocomments::Bool=false)
    status = Ref{Cint}(0)
    header = Ref{Ptr{UInt8}}()
    nkeys = Ref{Cint}(0)
    ccall((:ffhdr2str, libcfitsio), Cint,
          (Ptr{Void}, Cint, Ptr{Ptr{UInt8}}, Cint,
           Ptr{Ptr{UInt8}}, Ref{Cint}, Ref{Cint}),
          f.ptr, nocomments, &C_NULL, 0, header, nkeys, status)
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

for (a,b) in ((:fits_movabs_hdu,"ffmahd"),
              (:fits_movrel_hdu,"ffmrhd"))
    @eval begin
        function ($a)(f::FITSFile, hduNum::Integer)
            hdu_type = Ref{Cint}(0)
            status = Ref{Cint}(0)
            ccall(($b,libcfitsio), Cint,
                  (Ptr{Void}, Cint, Ref{Cint}, Ref{Cint}),
                  f.ptr, hduNum, hdu_type, status)
            fits_assert_ok(status[])
            hdu_int_to_type(hdu_type[])
        end
    end
end

function fits_movnam_hdu(f::FITSFile, extname::Compat.ASCIIString, extver::Integer=0,
                         hdu_type::Integer=-1)
    status = Ref{Cint}(0)
    ccall((:ffmnhd,libcfitsio), Cint,
          (Ptr{Void}, Cint, Ptr{UInt8}, Cint, Ref{Cint}),
          f.ptr, hdu_type, extname, extver, status)
    fits_assert_ok(status[])
end

function fits_get_hdu_num(f::FITSFile)
    hdunum = Ref{Cint}(0)
    ccall((:ffghdn,libcfitsio), Cint,
          (Ptr{Void}, Ref{Cint}),
          f.ptr, hdunum)
    hdunum[]
end

function fits_get_hdu_type(f::FITSFile)
    hdutype = Ref{Cint}(0)
    status = Ref{Cint}(0)
    ccall((:ffghdt, libcfitsio), Cint,
          (Ptr{Void}, Ref{Cint}, Ref{Cint}),
          f.ptr, hdutype, status)
    fits_assert_ok(status[])
    hdu_int_to_type(hdutype[])
end

# -----------------------------------------------------------------------------
# image HDU functions

for (a, b) in ((:fits_get_img_type,      "ffgidt"),
               (:fits_get_img_equivtype, "ffgiet"),
               (:fits_get_img_dim,       "ffgidm"))
    @eval function ($a)(f::FITSFile)
        result = Ref{Cint}(0)
        status = Ref{Cint}(0)
        ccall(($b, libcfitsio), Cint, (Ptr{Void}, Ref{Cint}, Ref{Cint}),
              f.ptr, result, status)
        fits_assert_ok(status[])
        result[]
    end
end

function fits_create_img{T, S<:Integer}(f::FITSFile, ::Type{T},
                                        naxes::Vector{S})
    status = Ref{Cint}(0)
    ccall((:ffcrimll, libcfitsio), Cint,
          (Ptr{Void}, Cint, Cint, Ptr{Int64}, Ref{Cint}),
          f.ptr, bitpix_from_type(T), length(naxes),
          Vector{Int64}(naxes), status)
    fits_assert_ok(status[])
end

function fits_write_pix{S<:Integer,T}(f::FITSFile, fpixel::Vector{S},
                                      nelements::Integer, data::Array{T})
    status = Ref{Cint}(0)
    ccall((:ffppxll, libcfitsio), Cint,
          (Ptr{Void}, Cint, Ptr{Int64}, Int64, Ptr{Void}, Ref{Cint}),
          f.ptr, cfitsio_typecode(T), Vector{Int64}(fpixel),
          nelements, data, status)
    fits_assert_ok(status[])
end

function fits_write_pix(f::FITSFile, data::Array)
    fits_write_pix(f, ones(Int64, length(size(data))), length(data), data)
end

function fits_read_pix{S<:Integer,T}(f::FITSFile, fpixel::Vector{S},
                                     nelements::Int, nullval::T,
                                     data::Array{T})
    anynull = Ref{Cint}(0)
    status = Ref{Cint}(0)
    ccall((:ffgpxvll, libcfitsio), Cint,
          (Ptr{Void}, Cint, Ptr{Int64}, Int64, Ptr{Void}, Ptr{Void},
           Ref{Cint}, Ref{Cint}),
          f.ptr, cfitsio_typecode(T), Vector{Int64}(fpixel),
          nelements, &nullval, data, anynull, status)
    fits_assert_ok(status[])
    anynull[]
end

function fits_read_pix{S<:Integer,T}(f::FITSFile, fpixel::Vector{S},
                                     nelements::Int, data::Array{T})
    anynull = Ref{Cint}(0)
    status = Ref{Cint}(0)
    ccall((:ffgpxvll, libcfitsio), Cint,
          (Ptr{Void}, Cint, Ptr{Int64}, Int64, Ptr{Void}, Ptr{Void},
           Ref{Cint}, Ref{Cint}),
          f.ptr, cfitsio_typecode(T), Vector{Int64}(fpixel),
          nelements, C_NULL, data, anynull, status)
    fits_assert_ok(status[])
    anynull[]
end

function fits_read_pix(f::FITSFile, data::Array)
    fits_read_pix(f, ones(Int64,length(size(data))), length(data), data)
end

function fits_read_subset{S1<:Integer,S2<:Integer,S3<:Integer,T}(
             f::FITSFile, fpixel::Vector{S1}, lpixel::Vector{S2},
             inc::Vector{S3}, data::Array{T})
    anynull = Ref{Cint}(0)
    status = Ref{Cint}(0)
    ccall((:ffgsv, libcfitsio), Cint,
          (Ptr{Void}, Cint, Ptr{Clong}, Ptr{Clong}, Ptr{Clong}, Ptr{Void},
           Ptr{Void}, Ref{Cint}, Ref{Cint}),
          f.ptr, cfitsio_typecode(T),
          Vector{Clong}(fpixel),
          Vector{Clong}(lpixel),
          Vector{Clong}(inc),
          C_NULL, data, anynull, status)
    fits_assert_ok(status[])
    anynull[]
end

function fits_copy_image_section(fin::FITSFile, fout::FITSFile,
                                 section::Compat.ASCIIString)
    status = Ref{Cint}(0)
    ccall((:fits_copy_image_section,libcfitsio), Cint,
          (Ptr{Void}, Ptr{Void}, Ptr{UInt8}, Ref{Cint}),
          fin.ptr, fout.ptr, section, status)
    fits_assert_ok(status[])
end


# -----------------------------------------------------------------------------
# ASCII/binary table HDU functions

# The three fields are: ttype, tform, tunit (CFITSIO's terminology)
const ColumnDef = Tuple{Compat.ASCIIString, Compat.ASCIIString, Compat.ASCIIString}

for (a,b) in ((:fits_create_binary_tbl, 2),
              (:fits_create_ascii_tbl,  1))
    @eval begin
        function ($a)(f::FITSFile, numrows::Integer,
                      coldefs::Array{ColumnDef}, extname::Compat.ASCIIString)

            # get length and convert coldefs to three arrays of Ptr{Uint8}
            ntype = length(coldefs)
            ttype = [pointer(x[1]) for x in coldefs]
            tform = [pointer(x[2]) for x in coldefs]
            tunit = [pointer(x[3]) for x in coldefs]
            status = Ref{Cint}(0)

            ccall(("ffcrtb", libcfitsio), Cint,
                  (Ptr{Void}, Cint, Int64, Cint,
                   Ptr{Ptr{UInt8}}, Ptr{Ptr{UInt8}},
                   Ptr{Ptr{UInt8}}, Ptr{UInt8}, Ref{Cint}),
                  f.ptr, $b, numrows, ntype,
                  ttype, tform, tunit, extname,
                  status)
            fits_assert_ok(status[])
        end
    end
end

for (a,b,T) in ((:fits_get_num_cols,  "ffgncl",  :Cint),
                (:fits_get_num_hdus,  "ffthdu",  :Cint),
                (:fits_get_rowsize,   "ffgrsz",  :Clong))
    @eval begin
        function ($a)(f::FITSFile)
            result = Ref{$T}(0)
            status = Ref{Cint}(0)
            ccall(($b,libcfitsio), Cint,
                  (Ptr{Void}, Ref{$T}, Ref{Cint}),
                  f.ptr, result, status)
            fits_assert_ok(status[])
            result[]
        end
    end
end

function fits_get_colnum(f::FITSFile, tmplt::Compat.ASCIIString)
    result = Ref{Cint}(0)
    status = Ref{Cint}(0)

    # Second argument is case-sensitivity of search: 0 = case-insensitive
    #                                                1 = case-sensitive
    ccall(("ffgcno", libcfitsio), Cint,
          (Ptr{Void}, Cint, Ptr{UInt8}, Ref{Cint}, Ref{Cint}),
          f.ptr, 0, tmplt, result, status)
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
@eval begin
    function fits_get_coltype(ff::FITSFile, colnum::Integer)
        typecode = Ref{Cint}(0)
        repcnt = Ref{$T}(0)
        width = Ref{$T}(0)
        status = Ref{Cint}(0)
        ccall(($ffgtcl,libcfitsio), Cint,
              (Ptr{Void}, Cint, Ref{Cint}, Ref{$T}, Ref{$T}, Ref{Cint}),
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
              (Ptr{Void}, Cint, Ref{Cint}, Ref{$T}, Ref{$T}, Ref{Cint}),
              ff.ptr, colnum, typecode, repcnt, width, status)
        fits_assert_ok(status[])
        return Int(typecode[]), Int(repcnt[]), Int(width[])
    end

    function fits_get_img_size(f::FITSFile)
        ndim = fits_get_img_dim(f)
        naxes = Vector{$T}(ndim)
        status = Ref{Cint}(0)
        ccall(($ffgisz, libcfitsio), Cint,
              (Ptr{Void}, Cint, Ptr{$T}, Ref{Cint}),
              f.ptr, ndim, naxes, status)
        fits_assert_ok(status[])
        naxes
    end

    function fits_get_num_rows(f::FITSFile)
        result = Ref{$T}(0)
        status = Ref{Cint}(0)
        ccall(($ffgnrw, libcfitsio), Cint,
              (Ptr{Void}, Ref{$T}, Ref{Cint}),
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
        naxes = Vector{$T}(99)  # 99 is the maximum allowed number of axes
        naxis = Ref{Cint}(0)
        status = Ref{Cint}(0)
        ccall(($ffgtdm,libcfitsio), Cint,
              (Ptr{Void}, Cint, Cint, Ref{Cint}, Ptr{$T}, Ref{Cint}),
              ff.ptr, colnum, length(naxes), naxis, naxes, status)
        fits_assert_ok(status[])
        return naxes[1:naxis[]]
    end

    function fits_write_tdim(ff::FITSFile, colnum::Integer,
                                 naxes::Array{$T})
        status = Ref{Cint}(0)
        ccall(($ffptdm, libcfitsio), Cint,
              (Ptr{Void}, Cint, Cint, Ptr{$T}, Ref{Cint}),
              ff.ptr, colnum, length(naxes), naxes, status)
        fits_assert_ok(status[])
    end

    function fits_read_descript(f::FITSFile, colnum::Integer, rownum::Integer)
        repeat = Ref{$T}(0)
        offset = Ref{$T}(0)
        status = Ref{Cint}(0)
        ccall(($ffgdes, libcfitsio), Cint,
              (Ptr{Void}, Cint, Int64, Ref{$T}, Ref{$T}, Ref{Cint}),
              f.ptr, colnum, rownum, repeat, offset, status)
        fits_assert_ok(status[])
        return Int(repeat[]), Int(offset[])
    end
end

function fits_read_col(f::FITSFile,
                       colnum::Integer,
                       firstrow::Integer,
                       firstelem::Integer,
                       data::Array{Compat.ASCIIString})

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
          (Ptr{Void}, Cint, Int64, Int64, Int64,
           Ptr{UInt8}, Ptr{Ptr{UInt8}}, Ref{Cint}, Ref{Cint}),
          f.ptr, colnum, firstrow, firstelem, length(data),
          " ", buffers, anynull, status)
    fits_assert_ok(status[])

    # Create strings out of the buffers, terminating at null characters.
    # Note that `Compat.ASCIIString(x)` does not copy the buffer x.
    for i in 1:length(data)
        zeropos = search(buffers[i], 0x00)
        data[i] = (zeropos >= 1) ? Compat.ASCIIString(buffers[i][1:(zeropos-1)]) :
                                   Compat.ASCIIString(buffers[i])
    end
end

function fits_read_col{T}(f::FITSFile,
                          colnum::Integer,
                          firstrow::Integer,
                          firstelem::Integer,
                          data::Array{T})
    anynull = Ref{Cint}(0)
    status = Ref{Cint}(0)
    ccall((:ffgcv,libcfitsio), Cint,
          (Ptr{Void}, Cint, Cint, Int64, Int64, Int64,
           Ptr{Void}, Ptr{Void}, Ref{Cint}, Ref{Cint}),
          f.ptr, cfitsio_typecode(T), colnum,
          firstrow, firstelem, length(data),
          C_NULL, data, anynull, status)
    fits_assert_ok(status[])
end

function fits_write_col(f::FITSFile,
                        colnum::Integer,
                        firstrow::Integer,
                        firstelem::Integer,
                        data::Array{Compat.ASCIIString})
    status = Ref{Cint}(0)
    ccall((:ffpcls, libcfitsio), Cint,
          (Ptr{Void}, Cint, Int64, Int64, Int64,
           Ptr{Ptr{UInt8}}, Ref{Cint}),
          f.ptr, colnum, firstrow, firstelem, length(data),
          data, status)
    fits_assert_ok(status[])
end

function fits_write_col{T}(f::FITSFile,
                           colnum::Integer,
                           firstrow::Integer,
                           firstelem::Integer,
                           data::Array{T})
    status = Ref{Cint}(0)
    ccall((:ffpcl, libcfitsio), Cint,
          (Ptr{Void}, Cint, Cint, Int64, Int64, Int64,
           Ptr{Void}, Ref{Cint}),
          f.ptr, cfitsio_typecode(T), colnum,
          firstrow, firstelem, length(data),
          data, status)
    fits_assert_ok(status[])
end

for (a,b) in ((:fits_insert_rows, "ffirow"),
              (:fits_delete_rows, "ffdrow"))
    @eval begin
        function ($a)(f::FITSFile, firstrow::Integer, nrows::Integer)
            status = Ref{Cint}(0)
            ccall(($b,libcfitsio), Cint,
                  (Ptr{Void}, Int64, Int64, Ref{Cint}),
                  f.ptr, firstrow, nrows, status)
            fits_assert_ok(status[])
        end
    end
end

# -----------------------------------------------------------------------------
# deprecations

import Base: @deprecate

# Deprecated in v0.7

@deprecate fits_get_num_rowsll fits_get_num_rows

# Deprecated in v0.6

@deprecate fits_get_col_repeat(f::FITSFile, colnum::Integer) fits_get_coltype(f, colnum)[2, 3]

@deprecate fits_read_col{T}(f::FITSFile, ::Type{T}, colnum::Integer, firstrow::Integer, firstelem::Integer, data::Array{T}) fits_read_col(f, colnum, firstrow, firstelem, data)

@deprecate fits_write_col{T}(f::FITSFile, ::Type{T}, colnum::Integer, firstrow::Integer, firstelem::Integer, data::Array{T}) fits_write_col(f, colnum, firstrow, firstelem, data)

end # module
