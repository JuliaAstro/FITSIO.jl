# API Reference

## File operations

```@docs
FITS
length(::FITS)
close(::FITS)
```

## Header operations

```@docs
read_header
FITSHeader
read_key
```

## Image operations

```@docs
write{T}(::FITS, ::Array{T})
read(::ImageHDU)
read(::ImageHDU, ::Union{Range{Int}, Int, Colon}...)
ndims(::ImageHDU)
size(::ImageHDU)
length(::ImageHDU)
copy_section
```

## Table operations

```@docs
write(::FITS, ::Dict{String})
write(::FITS, ::Vector{String}, ::Vector)
read(::TableHDU, ::String)
```

## Miscellaneous

```@docs
FITSIO.libcfitsio_version
```