# API Reference

## File operations

```@docs
FITS
length(::FITS)
close(::FITS)
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
```

## Image operations

```@docs
read(::ImageHDU)
read!
FITSIO.fitsread
write(::FITS, ::StridedArray{<:Real})
write(::ImageHDU, ::StridedArray{<:Real})
FITSIO.fitswrite
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
