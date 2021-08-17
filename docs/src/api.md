# API Reference

## File operations

```@docs
FITS
length
close
deleteat!
```

## Header operations

```@docs
read_key
write_key
read_header
FITSHeader
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
