# API Reference

## File operations

```@docs
FITS
length
flush
close
deleteat!
```

## HDU operations
```@docs
FITSIO.write_checksum
FITSIO.verify_checksum
```

## Header operations

```@docs
read_key
write_key
read_header
FITSHeader
FITSHeader(::AbstractVector{<:NamedTuple})
FITSHeader(::WCS.WCSTransform)
length(::FITSHeader)
haskey(::FITSHeader, ::String)
keys(::FITSHeader)
values(::FITSHeader)
get_comment
set_comment!
default_header
```

## Image operations

```@docs
read(::ImageHDU)
read!
FITSIO.fitsread
write(::FITS, ::StridedArray{<:Real})
write(::ImageHDU, ::StridedArray{<:Real})
FITSIO.fitswrite
eltype(::ImageHDU)
ndims(::ImageHDU)
size(::ImageHDU)
length(::ImageHDU)
copy_section
```

## Table operations

```@docs
FITSIO.colnames
write(::FITS, ::Dict{String})
write(::FITS, ::Vector{String}, ::Vector)
read(::TableHDU, ::String)
```
