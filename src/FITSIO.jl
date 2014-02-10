module FITSIO

export FITSFile,
       fits_clobber_file,
       fits_close_file,
       fits_create_file,
       fits_create_img,
       fits_delete_file,
       fits_delete_key,
       fits_delete_record,
       fits_file_mode,
       fits_file_name,
       fits_get_hdrspace,
       fits_get_hdu_num,
       fits_get_img_size,
       fits_get_num_cols,
       fits_get_num_hdus,
       fits_get_num_rows,
       fits_get_num_rowsll,
       fits_movabs_hdu,
       fits_movrel_hdu,
       fits_open_data,
       fits_open_file,
       fits_open_image,
       fits_open_table,
       fits_read_col,
       fits_read_keyn,
       fits_read_keyword,
       fits_read_pix,
       fits_read_record,
       fits_write_key,
       fits_write_pix,
       fits_write_record

export fitsread

import Base: close, show

type FITSFile
    ptr::Ptr{Void}
    status::Int32
end

function fits_get_errstatus(status::Int32)
    msg = Array(Uint8, 31)
    ccall((:ffgerr,:libcfitsio), Void, (Int32,Ptr{Uint8}), status, msg)
    bytestring(convert(Ptr{Uint8},msg))
end

function fits_assert_ok(status::Int32)
    if status != 0
        throw(fits_get_errstatus(status))
    end
end
fits_assert_ok(f::FITSFile) = fits_assert_ok(f.status)

# constants

_cfitsio_bitpix{T<:Integer}(::Type{T}) = int32(8*sizeof(T))
_cfitsio_bitpix{T<:FloatingPoint}(::Type{T}) = int32(-8*sizeof(T))

_cfitsio_datatype(::Type{Uint8})      = int32(11)
_cfitsio_datatype(::Type{Int8})       = int32(12)
_cfitsio_datatype(::Type{Bool})       = int32(14)
_cfitsio_datatype{T<:String}(::Type{T}) = int32(16)
_cfitsio_datatype(::Type{Uint16})     = int32(20)
_cfitsio_datatype(::Type{Int16})      = int32(21)
_cfitsio_datatype(::Type{Uint32})     = int32(30)
_cfitsio_datatype(::Type{Int32})      = int32(31)
_cfitsio_datatype(::Type{Float32})    = int32(42)
_cfitsio_datatype(::Type{Int64})      = int32(81)
_cfitsio_datatype(::Type{Float64})    = int32(82)
_cfitsio_datatype(::Type{Complex64})  = int32(83)
_cfitsio_datatype(::Type{Complex128}) = int32(163)

# General-purpose functions

for (a,b,T) in ((:fits_file_mode,     "ffflmd",  :Cint),
                (:fits_get_num_cols,  "ffgncl",  :Cint),
                (:fits_get_num_hdus,  "ffthdu",  :Cint),
                (:fits_get_num_rows,  "ffgnrw",  :Clong),
                (:fits_get_num_rowsll,"ffgnrwll",:Clonglong))
    @eval begin
        function ($a)(f::FITSFile)
            result = $T[0]
            ccall(($b,:libcfitsio), Int32,
                  (Ptr{Void}, Ptr{$T}, Ptr{Int32}),
                  f.ptr, result, &f.status)
            result[1]
        end
    end
end

# file access

function fits_create_file(filename::String)
    ptr = Array(Ptr{Void}, 1)
    status = Int32[0]
    ccall((:ffinit,:libcfitsio),
        Int32, (Ptr{Ptr{Void}},Ptr{Uint8},Ptr{Int32}),
        ptr, bytestring(filename), status)
    fits_assert_ok(status[1])
    FITSFile(ptr[1], status[1])
end

fits_clobber_file(filename::String) = fits_create_file("!"*filename)

for (a,b) in ((:fits_open_data, "ffdopn"),
              (:fits_open_file, "ffopen"),
              (:fits_open_image,"ffiopn"),
              (:fits_open_table,"fftopn"))
    @eval begin
        function ($a)(filename::String)
            ptr = Array(Ptr{Void}, 1)
            mode = int32(0) # readonly
            status = Int32[0]
            ccall(($b,:libcfitsio), Int32,
                  (Ptr{Ptr{Void}},Ptr{Uint8},Int32,Ptr{Int32}),
                  ptr, bytestring(filename), mode, status)
            fits_assert_ok(status[1])
            FITSFile(ptr[1], status[1])
        end
    end
end

for (a,b) in ((:fits_close_file, "ffclos"),
              (:fits_delete_file,"ffdelt"))
    @eval begin
        function ($a)(f::FITSFile)
            ccall(($b,:libcfitsio), Int32,
                  (Ptr{Void},Ptr{Int32}),
                  f.ptr, &f.status)
            fits_assert_ok(f)
        end
    end
end

close(f::FITSFile) = fits_close_file(f)

function fits_file_name(f::FITSFile)
    value = Array(Uint8, 1025)
    ccall((:ffflnm,:libcfitsio), Int32,
          (Ptr{Void},Ptr{Uint8},Ptr{Int32}),
          f.ptr, value, &f.status)
    fits_assert_ok(f)
    bytestring(convert(Ptr{Uint8}, value))
end

# header keywords

function fits_get_hdrspace(f::FITSFile)
    keysexist = Int32[0]
    morekeys = Int32[0]
    ccall((:ffghsp,:libcfitsio), Int32,
        (Ptr{Void},Ptr{Int32},Ptr{Int32},Ptr{Int32}),
        f.ptr, keysexist, morekeys, &f.status)
    (keysexist[1], morekeys[1])
end

function fits_read_keyword(f::FITSFile, keyname::String)
    value = Array(Uint8, 71)
    comment = Array(Uint8, 71)
    ccall((:ffgkey,:libcfitsio), Int32,
        (Ptr{Void},Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Ptr{Int32}),
        f.ptr, bytestring(keyname), value, comment, &f.status)
    fits_assert_ok(f)
    bytestring(convert(Ptr{Uint8},value)), bytestring(convert(Ptr{Uint8},comment))
end

function fits_read_record(f::FITSFile, keynum::Int)
    card = Array(Uint8, 81)
    ccall((:ffgrec,:libcfitsio), Int32,
        (Ptr{Void},Int32,Ptr{Uint8},Ptr{Int32}),
        f.ptr, keynum, card, &f.status)
    fits_assert_ok(f)
    bytestring(convert(Ptr{Uint8},card))
end

function fits_read_keyn(f::FITSFile, keynum::Int)
    keyname = Array(Uint8, 9)
    value = Array(Uint8, 71)
    comment = Array(Uint8, 71)
    ccall((:ffgkyn,:libcfitsio), Int32,
        (Ptr{Void},Int32,Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Ptr{Int32}),
        f.ptr, keynum, keyname, value, comment, &f.status)
    fits_assert_ok(f)
    bytestring(convert(Ptr{Uint8},keyname)),
    bytestring(convert(Ptr{Uint8},value)),
    bytestring(convert(Ptr{Uint8},comment))
end

function fits_write_key(f::FITSFile, keyname::String, value::Union(Number,String), comment::String)
    cvalue = isa(value,String) ?  bytestring(value) :
             isa(value,Bool) ? [int32(value)] : [value]
    ccall((:ffpky,:libcfitsio), Int32,
        (Ptr{Void},Int32,Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Ptr{Int32}),
        f.ptr, _cfitsio_datatype(typeof(value)), bytestring(keyname),
        cvalue, bytestring(comment), &f.status)
    fits_assert_ok(f)
end

function fits_write_record(f::FITSFile, card::String)
    ccall((:ffprec,:libcfitsio), Int32,
        (Ptr{Void},Ptr{Uint8},Ptr{Int32}),
        f.ptr, bytestring(card), &f.status)
    fits_assert_ok(f)
end

function fits_delete_record(f::FITSFile, keynum::Int)
    ccall((:ffdrec,:libcfitsio), Int32,
        (Ptr{Void},Int32,Ptr{Int32}),
        f.ptr, keynum, &f.status)
    fits_assert_ok(f)
end

function fits_delete_key(f::FITSFile, keyname::String)
    ccall((:ffdkey,:libcfitsio), Int32,
        (Ptr{Void},Ptr{Uint8},Ptr{Int32}),
        f.ptr, bytestring(keyname), &f.status)
    fits_assert_ok(f)
end

# HDU functions

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
            hdu_type = Int32[0]
            ccall(($b,:libcfitsio), Int32,
                  (Ptr{Void}, Int32, Ptr{Int32}, Ptr{Int32}),
                  f.ptr, hduNum, hdu_type, &f.status)
            fits_assert_ok(f)
            hdu_int_to_type(hdu_type[1])
        end
    end
end

function fits_get_hdu_num(f::FITSFile)
    hdunum = Int32[0]
    ccall((:ffghdn,:libcfitsio), Int32,
          (Ptr{Void},Ptr{Int32}),
          f.ptr, hdunum)
    hdunum[1]
end

# primary array or IMAGE extension

function fits_get_img_size(f::FITSFile)
    naxis = Int32[0]
    ccall((:ffgidm,:libcfitsio), Int32,
        (Ptr{Void},Ptr{Int32},Ptr{Int32}),
        f.ptr, naxis, &f.status)
    fits_assert_ok(f)
    naxes = zeros(Int, naxis[1])
    ccall((:ffgisz,:libcfitsio), Int32,
        (Ptr{Void},Int32,Ptr{Int},Ptr{Int32}),
        f.ptr, naxis[1], naxes, &f.status)
    fits_assert_ok(f)
    naxes
end

function fits_create_img(f::FITSFile, t::Type, naxes::Vector{Int})
    ccall((:ffcrim,:libcfitsio), Int32,
        (Ptr{Void},Int32,Int32,Ptr{Int},Ptr{Int32}),
        f.ptr, _cfitsio_bitpix(t), length(naxes), naxes, &f.status)
    fits_assert_ok(f)
end

function fits_write_pix{T}(f::FITSFile, fpixel::Vector{Int}, nelements::Int, data::Array{T})
    ccall((:ffppx,:libcfitsio), Int32,
        (Ptr{Void},Int32,Ptr{Int},Int,Ptr{Void},Ptr{Int32}),
        f.ptr, _cfitsio_datatype(T), fpixel, nelements, data, &f.status)
    fits_assert_ok(f)
end
fits_write_pix(f::FITSFile, data::Array) = fits_write_pix(f, ones(Int,length(size(data))), length(data), data)

function fits_read_pix{T}(f::FITSFile, fpixel::Vector{Int}, nelements::Int, nullval::T, data::Array{T})
    anynull = Int32[0]
    ccall((:ffgpxv,:libcfitsio), Int32,
        (Ptr{Void},Int32,Ptr{Int},Int,Ptr{Void},Ptr{Void},Ptr{Int32},Ptr{Int32}),
        f.ptr, _cfitsio_datatype(T), fpixel, nelements, &nullval, data, anynull, &f.status)
    fits_assert_ok(f)
    anynull[1]
end
function fits_read_pix{T}(f::FITSFile, fpixel::Vector{Int}, nelements::Int, data::Array{T})
    anynull = Int32[0]
    ccall((:ffgpxv,:libcfitsio), Int32,
        (Ptr{Void},Int32,Ptr{Int},Int,Ptr{Void},Ptr{Void},Ptr{Int32},Ptr{Int32}),
        f.ptr, _cfitsio_datatype(T), fpixel, nelements, C_NULL, data, anynull, &f.status)
    fits_assert_ok(f)
    anynull[1]
end
fits_read_pix(f::FITSFile, data::Array) = fits_read_pix(f, ones(Int,length(size(data))), length(data), data)

function fitsread(filename::String)
    f = fits_open_file(filename)
    s = fits_get_img_size(f)
    a = Array(Float64, s...)
    fits_read_pix(f, a)
    fits_close_file(f)
    a'
end

# ASCII/binary tables

function fits_read_col{T}(f::FITSFile,
                          ::Type{T},
                          colnum::Int,
                          firstrow::Int64,
                          firstelem::Int64,
                          nelements::Int64)

    result = zeros(T, nelements)
    anynull = Int32[0]
    nullvalue = T[0]

    ccall((:ffgcv,:libcfitsio), Int32,
          (Ptr{Void}, Int32, Int32, Int64, Int64, Int64,
           Ptr{T}, Ptr{T}, Ptr{Int32}, Ptr{Int32}),
          f.ptr, _cfitsio_datatype(T), convert(Int32, colnum),
          firstrow, firstelem, nelements,
          nullvalue, result, anynull, &f.status)

    return result

end

const mode_strs = [int32(0)=>"READONLY", int32(1)=>"READWRITE"]

function show(io::IO, f::FITSFile)
    print(io, "file: ", fits_file_name(f), "\n")
    print(io, "mode: ", mode_strs[fits_file_mode(f)], "\n")
    print(io, "extnum hdutype         hduname\n")

    current = fits_get_hdu_num(f)  # Mark the current HDU.

    for i = 1:fits_get_num_hdus(f)
        hdutype = fits_movabs_hdu(f, i)
        extname = fits_read_keyword(f, "EXTNAME")[1]
        @printf io "%-6d %-15s %s\n" i hdutype extname
    end
    fits_movabs_hdu(f, current)  # Return to the HDU we were on.
end

end # module
