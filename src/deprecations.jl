import Base: @deprecate

# Deprecated in v0.6

@deprecate fits_get_col_repeat(f::FITSFile, colnum::Integer) fits_get_coltype(f, colnum)[2, 3]

@deprecate fits_read_col{T}(f::FITSFile, ::Type{T}, colnum::Integer, firstrow::Integer, firstelem::Integer, data::Array{T}) fits_read_col(f, colnum, firstrow, firstelem, data)

@deprecate fits_write_col{T}(f::FITSFile, ::Type{T}, colnum::Integer, firstrow::Integer, firstelem::Integer, data::Array{T}) fits_write_col(f, colnum, firstrow, firstelem, data)
