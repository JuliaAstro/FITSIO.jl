using Compat

# functions to convert to C integer types (this could be in Julia Base module):
for (f, T) in ((:cint,  Cint),
               (:clong, Clong))
    @eval begin
        ($f)(x::$T) = x
        ($f)(x::Integer) = convert($T, x)
        ($f)(x::AbstractArray{$T}) = x
        function ($f){T<:Integer}(x::AbstractArray{T})
            y = similar(x, $T)
            i = 1
            for e in x
                y[i] = convert($T, e)
                i += 1
            end
            y
        end
    end
end

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
#       10  SBYTE_IMG       Int8     written as: BITPIX = 8, BSCALE = 1, BZERO = -128
#       20  USHORT_IMG      Uint16   written as: BITPIX = 16, BSCALE = 1, BZERO = 32768
#       40  LONG_IMG        Uint32   written as: BITPIX = 32, BSCALE = 1, BZERO = 2147483648
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

const bitpix_to_type = Dict{Cint, DataType}()
for (T, code) in ((Uint8,     8), # BYTE_IMG
                  (Int16,    16), # SHORT_IMG
                  (Int32,    32), # LONG_IMG
                  (Int64,    64), # LONGLONG_IMG
                  (Float32, -32), # FLOAT_IMG
                  (Float64, -64), # DOUBLE_IMG
                  (Int8,     10), # SBYTE_IMG
                  (Uint16,   20), # USHORT_IMG
                  (Uint32,   40)) # ULONG_IMG
    local value = cint(code)
    @eval begin
        bitpix_to_type[$value] = $T
        _cfitsio_bitpix(::Type{$T}) = $value
    end
end

for (T, code) in ((Uint8,       11),
                  (Int8,        12),
                  (Bool,        14),
                  (String,      16),
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
    local value = cint(code)
    @eval begin
        if $T == String
            _cfitsio_datatype{T<:$T}(::Type{T}) = $value
        else
            _cfitsio_datatype(::Type{$T}) = $value
        end
    end
end


type FITSFile
    ptr::Ptr{Void}

    function FITSFile(ptr::Ptr{Void})
        f = new(ptr)
        finalizer(f, fits_close_file)
        f
    end
end

function fits_assert_open(f::FITSFile)
    if f.ptr == C_NULL
        error("attempt to access closed FITS file")
    end
end

function fits_get_errstatus(status::Cint)
    msg = Array(Uint8, 31)
    ccall((:ffgerr,libcfitsio), Void, (Cint,Ptr{Uint8}), status, msg)
    bytestring(convert(Ptr{Uint8},msg))
end

function fits_assert_ok(status::Cint)
    if status != 0
        error(fits_get_errstatus(status))
    end
end

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


# General-purpose functions

for (a,b,T) in ((:fits_file_mode,     "ffflmd",  :Cint),
                (:fits_get_num_cols,  "ffgncl",  :Cint),
                (:fits_get_num_hdus,  "ffthdu",  :Cint),
                (:fits_get_num_rows,  "ffgnrw",  :Clong),
                (:fits_get_num_rowsll,"ffgnrwll",:Int64),
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

# file access

function fits_create_file(filename::String)
    ptr = Array(Ptr{Void}, 1)
    status = Cint[0]
    ccall((:ffinit,libcfitsio),
        Cint, (Ptr{Ptr{Void}},Ptr{Uint8},Ptr{Cint}),
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
    bytestring(convert(Ptr{Uint8}, value))
end

# header keywords

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

function fits_read_keyword(f::FITSFile, keyname::String)
    value = Array(Uint8, 71)
    comment = Array(Uint8, 71)
    status = Cint[0]
    ccall((:ffgkey,libcfitsio), Cint,
        (Ptr{Void},Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Ptr{Cint}),
        f.ptr, bytestring(keyname), value, comment, status)
    fits_assert_ok(status[1])
    bytestring(convert(Ptr{Uint8},value)),
    bytestring(convert(Ptr{Uint8},comment))
end

function fits_read_record(f::FITSFile, keynum::Integer)
    card = Array(Uint8, 81)
    status = Cint[0]
    ccall((:ffgrec,libcfitsio), Cint,
        (Ptr{Void},Cint,Ptr{Uint8},Ptr{Cint}),
        f.ptr, keynum, card, status)
    fits_assert_ok(status[1])
    bytestring(convert(Ptr{Uint8},card))
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
    bytestring(convert(Ptr{Uint8},keyname)),
    bytestring(convert(Ptr{Uint8},value)),
    bytestring(convert(Ptr{Uint8},comment))
end

function fits_write_key(f::FITSFile, keyname::String,
                        value::Union(FloatingPoint,String), comment::String)
    cvalue = isa(value,String) ?  bytestring(value) :
             isa(value,Bool) ? [cint(value)] : [value]
    status = Cint[0]
    ccall((:ffpky,libcfitsio), Cint,
        (Ptr{Void},Cint,Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Ptr{Cint}),
        f.ptr, _cfitsio_datatype(typeof(value)), bytestring(keyname),
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

function fits_write_record(f::FITSFile, card::String)
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

function fits_delete_key(f::FITSFile, keyname::String)
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

# HDU functions

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

function fits_movnam_hdu(f::FITSFile, extname::String, extver::Integer=0,
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

# primary array or IMAGE extension

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

function fits_create_img{S<:Integer}(f::FITSFile, t::Type, naxes::Vector{S})
    status = Cint[0]
    ccall((:ffcrimll, libcfitsio), Cint,
          (Ptr{Void}, Cint, Cint, Ptr{Int64}, Ptr{Cint}),
          f.ptr, _cfitsio_bitpix(t), length(naxes), int64(naxes), status)
    fits_assert_ok(status[1])
end

function fits_write_pix{S<:Integer,T}(f::FITSFile, fpixel::Vector{S},
                                      nelements::Integer, data::Array{T})
    status = Cint[0]
    ccall((:ffppxll, libcfitsio), Cint,
        (Ptr{Void}, Cint, Ptr{Int64}, Int64, Ptr{Void}, Ptr{Cint}),
        f.ptr, _cfitsio_datatype(T), int64(fpixel), nelements, data, status)
    fits_assert_ok(status[1])
end
fits_write_pix(f::FITSFile, data::Array) = fits_write_pix(f, ones(Int64, length(size(data))), length(data), data)

function fits_read_pix{S<:Integer,T}(f::FITSFile, fpixel::Vector{S}, nelements::Int, nullval::T, data::Array{T})
    anynull = Cint[0]
    status = Cint[0]
    ccall((:ffgpxv,libcfitsio), Cint,
          (Ptr{Void},Cint,Ptr{Clong},Int64,Ptr{Void},Ptr{Void},Ptr{Cint},Ptr{Cint}),
          f.ptr, _cfitsio_datatype(T), int64(fpixel), nelements, &nullval, data, anynull, status)
    fits_assert_ok(status[1])
    anynull[1]
end
function fits_read_pix{S<:Integer,T}(f::FITSFile, fpixel::Vector{S}, nelements::Int, data::Array{T})
    anynull = Cint[0]
    status = Cint[0]
    ccall((:ffgpxv,libcfitsio), Cint,
          (Ptr{Void},Cint,Ptr{Clong},Int64,Ptr{Void},Ptr{Void},Ptr{Cint},Ptr{Cint}),
          f.ptr, _cfitsio_datatype(T), int64(fpixel), nelements, C_NULL, data, anynull, status)
    fits_assert_ok(status[1])
    anynull[1]
end
fits_read_pix(f::FITSFile, data::Array) = fits_read_pix(f, ones(Int64,length(size(data))), length(data), data)

function fits_read_subset{S1<:Integer,S2<:Integer,S3<:Integer,T}(f::FITSFile,
                                                                 fpixel::Vector{S1},
                                                                 lpixel::Vector{S2},
                                                                 inc::Vector{S3},
                                                                 data::Array{T})
    anynull = Cint[0]
    status = Cint[0]
    ccall((:ffgsv,libcfitsio), Cint,
          (Ptr{Void},Cint,Ptr{Clong},Ptr{Clong},Ptr{Clong},Ptr{Void},Ptr{Void},Ptr{Cint},Ptr{Cint}),
          f.ptr, _cfitsio_datatype(T),
          clong(fpixel), clong(lpixel), clong(inc),
          C_NULL, data, anynull, status)
    fits_assert_ok(status[1])
    anynull[1]
end

function fits_copy_image_section(fin::FITSFile, fout::FITSFile,
                                 section::String)
    status = Cint[0]
    ccall((:fits_copy_image_section,libcfitsio), Cint,
          (Ptr{Void}, Ptr{Void}, Ptr{Uint8}, Ptr{Cint}),
          fin.ptr, fout.ptr, bytestring(section), status)
    fits_assert_ok(status[1])
end

# ASCII/binary tables

# The three fields are: ttype, tform, tunit (CFITSIO's terminology)
typealias ColumnDef (ASCIIString, ASCIIString, ASCIIString)

for (a,b) in ((:fits_create_binary_tbl, 2),
              (:fits_create_ascii_tbl,  1))
    @eval begin
        function ($a)(f::FITSFile, numrows::Integer,
                      coldefs::Array{ColumnDef}, extname::String)

            ntype = length(coldefs)
            ttype = map((x) -> pointer(x[1].data), coldefs)
            tform = map((x) -> pointer(x[2].data), coldefs)
            tunit = map((x) -> pointer(x[3].data), coldefs)
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

function fits_get_col_repeat(f::FITSFile, colnum::Integer)
    typecode = Cint[0]
    repeat = Int64[0]
    width = Int64[0]
    status = Cint[0]

    ccall((:ffgtclll, libcfitsio), Cint,
          (Ptr{Void}, Cint, Ptr{Cint}, Ptr{Int64}, Ptr{Int64}, Ptr{Cint}),
          f.ptr, colnum, typecode, repeat, width, status)

    (repeat[1], width[1])
end

function fits_read_col{T}(f::FITSFile,
                          ::Type{T},
                          colnum::Integer,
                          firstrow::Integer,
                          firstelem::Integer,
                          data::Array{T})

    anynull = Cint[0]
    status = Cint[0]
    nelements = length(data)

    if isa(T, Type{String}) || isa(T, Type{ASCIIString})

        # Make sure there is enough room for each string
        repcount, width = fits_get_col_repeat(f, colnum)
        for i in 1:length(data)
            # We need to call `repeat' N times in order to allocate N
            # strings
            data[i] = repeat(" ", repcount)
        end

        ccall((:ffgcvs, libcfitsio), Cint,
              (Ptr{Void}, Cint, Int64, Int64, Int64,
               Ptr{Uint8}, Ptr{Ptr{Uint8}}, Ptr{Cint}, Ptr{Cint}),
              f.ptr, colnum, firstrow, firstelem, nelements,
              "", data, anynull, status)

        # Truncate the strings to the first NULL character (if present)
        for idx in 1:length(data)
            zeropos = search(data[idx], '\0')
            if zeropos >= 1
                data[idx] = (data[idx])[1:(zeropos-1)]
            end
        end

    else

       ccall((:ffgcv,libcfitsio), Cint,
              (Ptr{Void}, Cint, Cint, Int64, Int64, Int64,
               Ptr{T}, Ptr{T}, Ptr{Cint}, Ptr{Cint}),
              f.ptr, _cfitsio_datatype(T), colnum,
              firstrow, firstelem, nelements,
              T[0], data, anynull, status)

    end

    fits_assert_ok(status[1])

end

function fits_write_col{T}(f::FITSFile,
                           ::Type{T},
                           colnum::Integer,
                           firstrow::Integer,
                           firstelem::Integer,
                           data::Array{T})

    status = Cint[0]
    nelements = length(data)

    if  isa(T, Type{String}) || isa(T, Type{ASCIIString})

        ccall((:ffpcls, libcfitsio), Cint,
              (Ptr{Void}, Cint, Int64, Int64, Int64,
               Ptr{Ptr{Uint8}}, Ptr{Cint}),
              f.ptr, colnum, firstrow, firstelem, nelements,
              data, status)

    else

        ccall((:ffpcl, libcfitsio), Cint,
              (Ptr{Void}, Cint, Cint, Int64, Int64, Int64,
               Ptr{T}, Ptr{Cint}),
              f.ptr, _cfitsio_datatype(T), colnum,
              firstrow, firstelem, nelements,
              data, status)

    end

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
