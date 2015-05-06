module FITSIO

using Compat

export FITS,
       HDU,
       ImageHDU,
       TableHDU,
       ASCIITableHDU,
       FITSHeader,
       read_key,
       read_header,
       get_comment,
       set_comment!,
       copy_section

import Base: getindex, setindex!, length, show, read, write, close, ndims,
             size, endof, haskey, keys, values

include("cfitsio.jl")  # Libcfitsio submodule
include("hdutypes.jl")
include("deprecations.jl")

end # module
