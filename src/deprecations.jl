# Deprecations in FITSIO module. Additional deprecations in FITSIO.Libcfitsio
# are in cfitsio.jl.

import Base: @deprecate

# Deprecated in v0.6
@deprecate readkey(hdu::HDU, key) read_key(hdu, key)
@deprecate readheader(hdu::HDU) read_header(hdu::HDU)
@deprecate getcomment(hdr::FITSHeader, key) get_comment(hdr, key)
@deprecate setcomment!(hdr::FITSHeader, key, comment) set_comment!(hdr, key, comment)
@deprecate getindex(hdu::ImageHDU, I...) read(hdu, I...)
