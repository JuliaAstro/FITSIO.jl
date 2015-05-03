import Base: @deprecate

# Deprecated in v0.6

export fits_get_col_repeat
@deprecate fits_get_col_repeat(f::FITSFile, colnum::Integer) fits_get_coltype(f, colnum)[2, 3]

@deprecate fits_read_col{T}(f::FITSFile, ::Type{T}, colnum::Integer, firstrow::Integer, firstelem::Integer, data::Array{T}) fits_read_col(f, colnum, firstrow, firstelem, data)

@deprecate fits_write_col{T}(f::FITSFile, ::Type{T}, colnum::Integer, firstrow::Integer, firstelem::Integer, data::Array{T}) fits_write_col(f, colnum, firstrow, firstelem, data)

export readkey, readheader, getcomment, setcomment!
@deprecate readkey(hdu::HDU, key) read_key(hdu, key)
@deprecate readheader(hdu::HDU) read_header(hdu::HDU)
@deprecate getcomment(hdr::FITSHeader, key) get_comment(hdr, key)
@deprecate setcomment!(hdr::FITSHeader, key, comment) set_comment!(hdr, key, comment)
@deprecate getindex(hdu::ImageHDU, I...) read(hdu, I...)
