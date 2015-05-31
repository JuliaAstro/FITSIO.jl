# cfitsio.jl - C-style interface to CFITSIO functions:
#
# - Function names closely mirror the C interface (e.g., `fits_open_file()`).
# - Functions operate on `FITSFile`, a thin wrapper for `fitsfile` C struct
#   (`FITSFile` has concept of "current HDU", as in CFITSIO).
# - Note that the wrapper functions *do* check the return status from CFITSIO
#   and throw an error with the appropriate message.
#
#
# Syntax compatibility notes:
#
# - `convert(Vector{T}, x)` can be changed to `Vector{T}(x)` once
#   v0.3 is no longer supported, or can be changed to
#   `@compat Vector{T}(x)` once syntax support is in Compat.
# - `convert(Cint, x)` can be changed to `Cint(x)` once v0.3 is not supported
#   or `@compat Cint(x)` once syntax support is in Compat.
# - Once v0.3 is no longer supported, the type `Nothing` should be changed to
#   `Void` (this is the type of `nothing`).
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
#       16  TSTRING         ASCIIString
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

module Libcfitsio

using Compat

export FITSFile,
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
       fits_read_col,
       fits_read_descript,
       fits_read_keyn,
       fits_read_keyword,
       fits_read_pix,
       fits_read_record,
       fits_read_subset,
       fits_update_key,
       fits_write_col,
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
    local value = convert(Cint, code)
    @eval begin
        TYPE_FROM_BITPIX[$value] = $T
        bitpix_from_type(::Type{$T}) = $value
    end
end

# TODO: elimiate duplicates (Clong vs Int64 or Clong vs Cint) ?
for (T, code) in ((UInt8,       11),
                  (Int8,        12),
                  (Bool,        14),
                  (ASCIIString, 16),
                  (Cushort,     20),
                  (Cshort,      21),
                  (Cuint,       30),
                  (Cint,        31),
                  (Culong,      40),
                  (Clong,       41),
                  (Float32,     42),
                  (Int64,       81),
                  (Float64,     82),
                  (Complex64,   83),
                  (Complex128, 163))
    @eval cfitsio_typecode(::Type{$T}) = convert(Cint, $code)
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

# -----------------------------------------------------------------------------
# error messaging

function fits_assert_open(f::FITSFile)
    if f.ptr == C_NULL
        error("attempt to access closed FITS file")
    end
end

function fits_get_errstatus(status::Cint)
    msg = Array(Uint8, 31)
    ccall((:ffgerr,libcfitsio), Void, (Cint,Ptr{Uint8}), status, msg)
    bytestring(pointer(msg))
end

function fits_assert_ok(status::Cint)
    if status != 0
        error(fits_get_errstatus(status))
    end
end

fits_get_version() = ccall((:ffvers, libcfitsio), Cfloat, (Ptr{Cfloat},), &0.)

# -----------------------------------------------------------------------------
# file access & info functions

function fits_create_file(filename::String)
    ptr = Array(Ptr{Void}, 1)
    status = Cint[0]
    ccall((:ffinit,libcfitsio), Cint, (Ptr{Ptr{Void}},Ptr{Uint8},Ptr{Cint}),
          ptr, bytestring(filename), status)
    fits_assert_ok(status[1])
    FITSFile(ptr[1])
end

fits_clobber_file(filename::String) = fits_create_file("!"*filename)

for (a,b) in ((:fits_open_data, "ffdopn"),
              (:fits_open_file, "ffopen"),
              (:fits_open_image,"ffiopn"),
              (:fits_open_table,"fftopn"))
    @eval begin
        function ($a)(filename::String, mode::Integer=0)
            ptr = Array(Ptr{Void}, 1)
            status = Cint[0]
            ccall(($b,libcfitsio), Cint,
                  (Ptr{Ptr{Void}},Ptr{Uint8},Cint,Ptr{Cint}),
                  ptr, bytestring(filename), mode, status)
            fits_assert_ok(status[1])
            FITSFile(ptr[1])
        end
    end
end

for (a,b) in ((:fits_close_file, "ffclos"),
              (:fits_delete_file,"ffdelt"))
    @eval begin
        function ($a)(f::FITSFile)

            # fits_close_file() is called during garbage collection, but file
            # may already be closed by user, so we need to check if it is open.
            if f.ptr != C_NULL
                status = Cint[0]
                ccall(($b,libcfitsio), Cint,
                      (Ptr{Void},Ptr{Cint}),
                      f.ptr, status)
                fits_assert_ok(status[1])
                f.ptr = C_NULL
            end
        end
    end
end

close(f::FITSFile) = fits_close_file(f)

function fits_file_name(f::FITSFile)
    value = Array(Uint8, 1025)
    status = Cint[0]
    ccall((:ffflnm,libcfitsio), Cint,
          (Ptr{Void},Ptr{Uint8},Ptr{Cint}),
          f.ptr, value, status)
    fits_assert_ok(status[1])
    bytestring(pointer(value))
end

function fits_file_mode(f::FITSFile)
    result = Cint[0]
    status = Cint[0]
    ccall(("ffflmd", libcfitsio), Cint, (Ptr{Void}, Ptr{Cint}, Ptr{Cint}),
          f.ptr, result, status)
    fits_assert_ok(status[1])
    result[1]
end


# -----------------------------------------------------------------------------
# header access functions

function fits_get_hdrspace(f::FITSFile)
    keysexist = Cint[0]
    morekeys = Cint[0]
    status = Cint[0]
    ccall((:ffghsp,libcfitsio), Cint,
        (Ptr{Void},Ptr{Cint},Ptr{Cint},Ptr{Cint}),
        f.ptr, keysexist, morekeys, status)
    fits_assert_ok(status[1])
    (keysexist[1], morekeys[1])
end

function fits_read_keyword(f::FITSFile, keyname::ASCIIString)
    value = Array(Uint8, 71)
    comment = Array(Uint8, 71)
    status = Cint[0]
    ccall((:ffgkey,libcfitsio), Cint,
        (Ptr{Void},Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Ptr{Cint}),
        f.ptr, bytestring(keyname), value, comment, status)
    fits_assert_ok(status[1])
    bytestring(pointer(value)), bytestring(pointer(comment))
end

function fits_read_record(f::FITSFile, keynum::Integer)
    card = Array(Uint8, 81)
    status = Cint[0]
    ccall((:ffgrec,libcfitsio), Cint,
        (Ptr{Void},Cint,Ptr{Uint8},Ptr{Cint}),
        f.ptr, keynum, card, status)
    fits_assert_ok(status[1])
    bytestring(pointer(card))
end

function fits_read_keyn(f::FITSFile, keynum::Integer)
    keyname = Array(Uint8, 9)
    value = Array(Uint8, 71)
    comment = Array(Uint8, 71)
    status = Cint[0]
    ccall((:ffgkyn,libcfitsio), Cint,
        (Ptr{Void},Cint,Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Ptr{Cint}),
        f.ptr, keynum, keyname, value, comment, status)
    fits_assert_ok(status[1])
    (bytestring(pointer(keyname)), bytestring(pointer(value)),
     bytestring(pointer(comment)))
end

function fits_write_key(f::FITSFile, keyname::ASCIIString,
                        value::Union(FloatingPoint,ASCIIString),
                        comment::ASCIIString)
    cvalue = isa(value,ASCIIString) ?  bytestring(value) :
             isa(value,Bool) ? Cint[value] : [value]
    status = Cint[0]
    ccall((:ffpky,libcfitsio), Cint,
        (Ptr{Void},Cint,Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Ptr{Cint}),
        f.ptr, cfitsio_typecode(typeof(value)), bytestring(keyname),
        cvalue, bytestring(comment), status)
    fits_assert_ok(status[1])
end

function fits_write_comment(f::FITSFile, comment::ASCIIString)
    status = Cint[0]
    ccall((:ffpcom, libcfitsio), Cint, (Ptr{Void}, Ptr{Uint8}, Ptr{Cint}),
          f.ptr, bytestring(comment), status)
    fits_assert_ok(status[1])
end

function fits_write_history(f::FITSFile, history::ASCIIString)
    status = Cint[0]
    ccall((:ffphis, libcfitsio), Cint, (Ptr{Void}, Ptr{Uint8}, Ptr{Cint}),
          f.ptr, bytestring(history), status)
    fits_assert_ok(status[1])
end

# update key: if already present, update it, otherwise add it.
for (a,T,S) in (("ffukys", :ASCIIString, :(Ptr{Uint8})),
                ("ffukyl", :Bool,        :Cint),
                ("ffukyj", :Integer,     :Int64))
    @eval begin
        function fits_update_key(f::FITSFile, key::ASCIIString, value::$T,
                                 comment::Union(ASCIIString, Ptr{Void})=C_NULL)
            status = Cint[0]
            ccall(($a, libcfitsio), Cint,
                  (Ptr{Void}, Ptr{Uint8}, $S, Ptr{Uint8}, Ptr{Cint}),
                  f.ptr, key, value, comment, status)
            fits_assert_ok(status[1])
        end
    end
end

function fits_update_key(f::FITSFile, key::ASCIIString, value::FloatingPoint,
                         comment::Union(ASCIIString, Ptr{Void})=C_NULL)
    status = Cint[0]
    ccall(("ffukyd", libcfitsio), Cint,
          (Ptr{Void}, Ptr{Uint8}, Cdouble, Cint, Ptr{Uint8}, Ptr{Cint}),
          f.ptr, key, value, -15, comment, status)
    fits_assert_ok(status[1])
end

function fits_update_key(f::FITSFile, key::ASCIIString, value::Nothing,
                         comment::Union(ASCIIString, Ptr{Void})=C_NULL)
    status = Cint[0]
    ccall(("ffukyu", libcfitsio), Cint,
          (Ptr{Void}, Ptr{Uint8}, Ptr{Uint8}, Ptr{Cint}),
          f.ptr, key, comment, status)
    fits_assert_ok(status[1])
end

function fits_write_record(f::FITSFile, card::ASCIIString)
    status = Cint[0]
    ccall((:ffprec,libcfitsio), Cint,
        (Ptr{Void},Ptr{Uint8},Ptr{Cint}),
        f.ptr, bytestring(card), status)
    fits_assert_ok(status[1])
end

function fits_delete_record(f::FITSFile, keynum::Integer)
    status = Cint[0]
    ccall((:ffdrec,libcfitsio), Cint,
        (Ptr{Void},Cint,Ptr{Cint}),
        f.ptr, keynum, status)
    fits_assert_ok(status[1])
end

function fits_delete_key(f::FITSFile, keyname::ASCIIString)
    status = Cint[0]
    ccall((:ffdkey,libcfitsio), Cint,
        (Ptr{Void},Ptr{Uint8},Ptr{Cint}),
        f.ptr, bytestring(keyname), status)
    fits_assert_ok(status[1])
end

function fits_hdr2str(f::FITSFile, nocomments::Bool=false)
    status = Cint[0]
    header = Array(Ptr{Uint8}, 1)
    nkeys = Cint[0]
    ccall((:ffhdr2str, libcfitsio), Cint,
          (Ptr{Void}, Cint, Ptr{Ptr{Uint8}}, Cint,
           Ptr{Ptr{Uint8}}, Ptr{Cint}, Ptr{Cint}),
          f.ptr, nocomments, &C_NULL, 0, header, nkeys, status)
    result = bytestring(header[1])

    # free header pointer allocated by cfitsio (result is a copy)
    ccall((:fffree, libcfitsio), Ptr{Cint},
          (Ptr{Uint8}, Ptr{Cint}),
          header[1], status)
    fits_assert_ok(status[1])
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
            hdu_type = Cint[0]
            status = Cint[0]
            ccall(($b,libcfitsio), Cint,
                  (Ptr{Void}, Cint, Ptr{Cint}, Ptr{Cint}),
                  f.ptr, hduNum, hdu_type, status)
            fits_assert_ok(status[1])
            hdu_int_to_type(hdu_type[1])
        end
    end
end

function fits_movnam_hdu(f::FITSFile, extname::ASCIIString, extver::Integer=0,
                         hdu_type::Integer=-1)
    status = Cint[0]
    ccall((:ffmnhd,libcfitsio), Cint,
          (Ptr{Void}, Cint, Ptr{Uint8}, Cint, Ptr{Cint}),
          f.ptr, hdu_type, bytestring(extname), extver, status)
    fits_assert_ok(status[1])
end

function fits_get_hdu_num(f::FITSFile)
    hdunum = Cint[0]
    ccall((:ffghdn,libcfitsio), Cint,
          (Ptr{Void}, Ptr{Cint}),
          f.ptr, hdunum)
    hdunum[1]
end

function fits_get_hdu_type(f::FITSFile)
    hdutype = Cint[0]
    status = Cint[0]
    ccall((:ffghdt, libcfitsio), Cint,
          (Ptr{Void}, Ptr{Cint}, Ptr{Cint}),
          f.ptr, hdutype, status)
    fits_assert_ok(status[1])
    hdu_int_to_type(hdutype[1])
end

# -----------------------------------------------------------------------------
# image HDU functions

for (a, b) in ((:fits_get_img_type,      "ffgidt"),
               (:fits_get_img_equivtype, "ffgiet"),
               (:fits_get_img_dim,       "ffgidm"))
    @eval function ($a)(f::FITSFile)
        result = Cint[0]
        status = Cint[0]
        ccall(($b, libcfitsio), Cint, (Ptr{Void}, Ptr{Cint}, Ptr{Cint}),
              f.ptr, result, status)
        fits_assert_ok(status[1])
        result[1]
    end
end

function fits_get_img_size(f::FITSFile)
    ndim = fits_get_img_dim(f)
    naxes = Array(Clong, ndim)
    status = Cint[0]
    ccall((:ffgisz, libcfitsio), Cint,
        (Ptr{Void}, Cint, Ptr{Clong}, Ptr{Cint}),
        f.ptr, ndim, naxes, status)
    fits_assert_ok(status[1])
    naxes
end

function fits_create_img{T, S<:Integer}(f::FITSFile, ::Type{T},
                                        naxes::Vector{S})
    status = Cint[0]
    ccall((:ffcrimll, libcfitsio), Cint,
          (Ptr{Void}, Cint, Cint, Ptr{Int64}, Ptr{Cint}),
          f.ptr, bitpix_from_type(T), length(naxes),
          convert(Vector{Int64}, naxes), status)
    fits_assert_ok(status[1])
end

function fits_write_pix{S<:Integer,T}(f::FITSFile, fpixel::Vector{S},
                                      nelements::Integer, data::Array{T})
    status = Cint[0]
    ccall((:ffppxll, libcfitsio), Cint,
          (Ptr{Void}, Cint, Ptr{Int64}, Int64, Ptr{Void}, Ptr{Cint}),
          f.ptr, cfitsio_typecode(T), convert(Vector{Int64}, fpixel),
          nelements, data, status)
    fits_assert_ok(status[1])
end

function fits_write_pix(f::FITSFile, data::Array)
    fits_write_pix(f, ones(Int64, length(size(data))), length(data), data)
end

function fits_read_pix{S<:Integer,T}(f::FITSFile, fpixel::Vector{S},
                                     nelements::Int, nullval::T,
                                     data::Array{T})
    anynull = Cint[0]
    status = Cint[0]
    ccall((:ffgpxv, libcfitsio), Cint,
          (Ptr{Void}, Cint, Ptr{Clong}, Int64, Ptr{Void}, Ptr{Void},
           Ptr{Cint}, Ptr{Cint}),
          f.ptr, cfitsio_typecode(T), convert(Vector{Int64}, fpixel),
          nelements, &nullval, data, anynull, status)
    fits_assert_ok(status[1])
    anynull[1]
end

function fits_read_pix{S<:Integer,T}(f::FITSFile, fpixel::Vector{S},
                                     nelements::Int, data::Array{T})
    anynull = Cint[0]
    status = Cint[0]
    ccall((:ffgpxv, libcfitsio), Cint,
          (Ptr{Void}, Cint, Ptr{Clong}, Int64, Ptr{Void}, Ptr{Void},
           Ptr{Cint}, Ptr{Cint}),
          f.ptr, cfitsio_typecode(T), convert(Vector{Int64}, fpixel),
          nelements, C_NULL, data, anynull, status)
    fits_assert_ok(status[1])
    anynull[1]
end

function fits_read_pix(f::FITSFile, data::Array)
    fits_read_pix(f, ones(Int64,length(size(data))), length(data), data)
end

function fits_read_subset{S1<:Integer,S2<:Integer,S3<:Integer,T}(
             f::FITSFile, fpixel::Vector{S1}, lpixel::Vector{S2},
             inc::Vector{S3}, data::Array{T})
    anynull = Cint[0]
    status = Cint[0]
    ccall((:ffgsv, libcfitsio), Cint,
          (Ptr{Void}, Cint, Ptr{Clong}, Ptr{Clong}, Ptr{Clong}, Ptr{Void},
           Ptr{Void}, Ptr{Cint}, Ptr{Cint}),
          f.ptr, cfitsio_typecode(T),
          convert(Vector{Clong}, fpixel),
          convert(Vector{Clong}, lpixel),
          convert(Vector{Clong}, inc),
          C_NULL, data, anynull, status)
    fits_assert_ok(status[1])
    anynull[1]
end

function fits_copy_image_section(fin::FITSFile, fout::FITSFile,
                                 section::ASCIIString)
    status = Cint[0]
    ccall((:fits_copy_image_section,libcfitsio), Cint,
          (Ptr{Void}, Ptr{Void}, Ptr{Uint8}, Ptr{Cint}),
          fin.ptr, fout.ptr, bytestring(section), status)
    fits_assert_ok(status[1])
end


# -----------------------------------------------------------------------------
# ASCII/binary table HDU functions

# The three fields are: ttype, tform, tunit (CFITSIO's terminology)
typealias ColumnDef @compat Tuple{ASCIIString, ASCIIString, ASCIIString}

for (a,b) in ((:fits_create_binary_tbl, 2),
              (:fits_create_ascii_tbl,  1))
    @eval begin
        function ($a)(f::FITSFile, numrows::Integer,
                      coldefs::Array{ColumnDef}, extname::ASCIIString)

            # get length and convert coldefs to three arrays of Ptr{Uint8}
            ntype = length(coldefs)
            ttype = [pointer(x[1]) for x in coldefs]
            tform = [pointer(x[2]) for x in coldefs]
            tunit = [pointer(x[3]) for x in coldefs]
            status = Cint[0]

            ccall(("ffcrtb", libcfitsio), Cint,
                  (Ptr{Void}, Cint, Int64, Cint,
                   Ptr{Ptr{Uint8}}, Ptr{Ptr{Uint8}},
                   Ptr{Ptr{Uint8}}, Ptr{Uint8}, Ptr{Cint}),
                  f.ptr, $b, numrows, ntype,
                  ttype, tform, tunit, bytestring(extname),
                  status)
            fits_assert_ok(status[1])
        end
    end
end

for (a,b,T) in ((:fits_get_num_cols,  "ffgncl",  :Cint),
                (:fits_get_num_hdus,  "ffthdu",  :Cint),
                (:fits_get_rowsize,   "ffgrsz",  :Clong))
    @eval begin
        function ($a)(f::FITSFile)
            result = $T[0]
            status = Cint[0]
            ccall(($b,libcfitsio), Cint,
                  (Ptr{Void}, Ptr{$T}, Ptr{Cint}),
                  f.ptr, result, status)
            fits_assert_ok(status[1])
            result[1]
        end
    end
end

function fits_get_colnum(f::FITSFile, tmplt::ASCIIString)
    result = Cint[0]
    status = Cint[0]

    # Second argument is case-sensitivity of search: 0 = case-insensitive
    #                                                1 = case-sensitive
    ccall(("ffgcno", libcfitsio), Cint,
          (Ptr{Void}, Cint, Ptr{Uint8}, Ptr{Cint}, Ptr{Cint}),
          f.ptr, 0, tmplt, result, status)
    fits_assert_ok(status[1])
    return result[1]
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
else
    T = Int64
    ffgtdm = "ffgtdmll"
    ffgnrw = "ffgnrwll"
    ffptdm = "ffptdmll"
    ffgtcl = "ffgtclll"
    ffeqty = "ffeqtyll"
    ffgdes = "ffgdesll"
end
@eval begin
    function fits_get_coltype(ff::FITSFile, colnum::Integer)
        typecode = Cint[0]
        repcnt = $T[0]
        width = $T[0]
        status = Cint[0]
        ccall(($ffgtcl,libcfitsio), Cint,
              (Ptr{Void}, Cint, Ptr{Cint}, Ptr{$T}, Ptr{$T}, Ptr{Cint}),
              ff.ptr, colnum, typecode, repcnt, width, status)
        fits_assert_ok(status[1])
        return @compat Int(typecode[1]), Int(repcnt[1]), Int(width[1])
    end

    function fits_get_eqcoltype(ff::FITSFile, colnum::Integer)
        typecode = Cint[0]
        repcnt = $T[0]
        width = $T[0]
        status = Cint[0]
        ccall(($ffeqty,libcfitsio), Cint,
              (Ptr{Void}, Cint, Ptr{Cint}, Ptr{$T}, Ptr{$T}, Ptr{Cint}),
              ff.ptr, colnum, typecode, repcnt, width, status)
        fits_assert_ok(status[1])
        return @compat Int(typecode[1]), Int(repcnt[1]), Int(width[1])
    end

    function fits_get_num_rows(f::FITSFile)
        result = $T[0]
        status = Cint[0]
        ccall(($ffgnrw, libcfitsio), Cint,
              (Ptr{Void}, Ptr{$T}, Ptr{Cint}),
              f.ptr, result, status)
        fits_assert_ok(status[1])
        return @compat Int(result[1])
    end

    # `fits_read_tdim` returns the dimensions of a table column in a
    # binary table. Normally this information is given by the TDIMn
    # keyword, but if this keyword is not present then this routine
    # returns `[r]` with `r` equals to the repeat count in the TFORM
    # keyword.
    function fits_read_tdim(ff::FITSFile, colnum::Integer)
        naxes = Array($T, 99)  # 99 is the maximum allowed number of axes
        naxis = Cint[0]
        status = Cint[0]
        ccall(($ffgtdm,libcfitsio), Cint,
              (Ptr{Void}, Cint, Cint, Ptr{Cint}, Ptr{$T}, Ptr{Cint}),
              ff.ptr, colnum, length(naxes), naxis, naxes, status)
        fits_assert_ok(status[1])
        return naxes[1:naxis[1]]
    end

    function fits_write_tdim(ff::FITSFile, colnum::Integer,
                                 naxes::Array{$T})
        status = Cint[0]
        ccall(($ffptdm, libcfitsio), Cint,
              (Ptr{Void}, Cint, Cint, Ptr{$T}, Ptr{Cint}),
              ff.ptr, colnum, length(naxes), naxes, status)
        fits_assert_ok(status[1])
    end

    function fits_read_descript(f::FITSFile, colnum::Integer, rownum::Integer)
        repeat = $T[0]
        offset = $T[0]
        status = Cint[0]
        ccall(($ffgdes, libcfitsio), Cint,
              (Ptr{Void}, Cint, Int64, Ptr{$T}, Ptr{$T}, Ptr{Cint}),
              f.ptr, colnum, rownum, repeat, offset, status)
        fits_assert_ok(status[1])
        return @compat Int(repeat[1]), Int(offset[1])
    end
end

function fits_read_col(f::FITSFile,
                       colnum::Integer,
                       firstrow::Integer,
                       firstelem::Integer,
                       data::Array{ASCIIString})

    # get width: number of characters in each string
    typecode, repcount, width = fits_get_eqcoltype(f, colnum)

    # ensure that data are strings, otherwise cfitsio will try to write
    # formatted strings, which have widths given by fits_get_col_display_width
    # not by the repeat value from fits_get_coltype.
    abs(typecode) == 16 || error("not a string column")

    # create an array of character buffers of the correct width
    buffers = [Array(Uint8, width) for i in 1:length(data)]

    # Call the CFITSIO function
    anynull = Cint[0]
    status = Cint[0]
    ccall((:ffgcvs, libcfitsio), Cint,
          (Ptr{Void}, Cint, Int64, Int64, Int64,
           Ptr{Uint8}, Ptr{Ptr{Uint8}}, Ptr{Cint}, Ptr{Cint}),
          f.ptr, colnum, firstrow, firstelem, length(data),
          " ", buffers, anynull, status)
    fits_assert_ok(status[1])

    # Create strings out of the buffers, terminating at null characters.
    # Note that `ASCIIString(x)` does not copy the buffer x.
    for i in 1:length(data)
        zeropos = search(buffers[i], 0x00)
        data[i] = (zeropos >= 1) ? ASCIIString(buffers[i][1:(zeropos-1)]) :
                                   ASCIIString(buffers[i])
    end
end

function fits_read_col{T}(f::FITSFile,
                          colnum::Integer,
                          firstrow::Integer,
                          firstelem::Integer,
                          data::Array{T})
    anynull = Cint[0]
    status = Cint[0]
    ccall((:ffgcv,libcfitsio), Cint,
          (Ptr{Void}, Cint, Cint, Int64, Int64, Int64,
           Ptr{Void}, Ptr{Void}, Ptr{Cint}, Ptr{Cint}),
          f.ptr, cfitsio_typecode(T), colnum,
          firstrow, firstelem, length(data),
          C_NULL, data, anynull, status)
    fits_assert_ok(status[1])
end

function fits_write_col(f::FITSFile,
                        colnum::Integer,
                        firstrow::Integer,
                        firstelem::Integer,
                        data::Array{ASCIIString})
    status = Cint[0]
    ccall((:ffpcls, libcfitsio), Cint,
          (Ptr{Void}, Cint, Int64, Int64, Int64,
           Ptr{Ptr{Uint8}}, Ptr{Cint}),
          f.ptr, colnum, firstrow, firstelem, length(data),
          data, status)
    fits_assert_ok(status[1])
end

function fits_write_col{T}(f::FITSFile,
                           colnum::Integer,
                           firstrow::Integer,
                           firstelem::Integer,
                           data::Array{T})
    status = Cint[0]
    ccall((:ffpcl, libcfitsio), Cint,
          (Ptr{Void}, Cint, Cint, Int64, Int64, Int64,
           Ptr{Void}, Ptr{Cint}),
          f.ptr, cfitsio_typecode(T), colnum,
          firstrow, firstelem, length(data),
          data, status)
    fits_assert_ok(status[1])
end

for (a,b) in ((:fits_insert_rows, "ffirow"),
              (:fits_delete_rows, "ffdrow"))
    @eval begin
        function ($a)(f::FITSFile, firstrow::Integer, nrows::Integer)
            status = Cint[0]
            ccall(($b,libcfitsio), Cint,
                  (Ptr{Void}, Int64, Int64, Ptr{Cint}),
                  f.ptr, firstrow, nrows, status)
            fits_assert_ok(status[1])
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
