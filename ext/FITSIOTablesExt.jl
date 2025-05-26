module FITSIOTablesExt

import Tables
using FITSIO: ASCIITableHDU, TableHDU, fits_get_col_info

#Tables.jl integration

const EitherTableHDU = Union{TableHDU, ASCIITableHDU}
Tables.istable(::Type{<:EitherTableHDU}) = true
Tables.columnaccess(::Type{<:EitherTableHDU}) = true
Tables.columns(t::EitherTableHDU) = t

function Tables.columnnames(t::EitherTableHDU)
    cns = FITSIO.colnames(t)
    #filter out variable-length cols
    isvar = [last(fits_get_col_info(t.fitsfile, i)) for i in 1:length(cns)]
    Symbol.(cns[.! isvar])
end
function Tables.getcolumn(t::EitherTableHDU, s::Symbol)
    col = FITSIO.read(t, String(s))
    if fits_get_col_info(t.fitsfile, findfirst(FITSIO.colnames(t) .== String(s)))[3]
        error("variable-length columns not supported in Tables.jl interface")
    end
    #reshape multidimensional array into array of (possibly multidimensional) arrays if necessary
    if (dim = length(size(col))) == 1
        col
    else
        [col[vcat(repeat([:], dim-1), i)...] for i in 1:size(col, dim)]
    end
end
#this uses Tables.columnnames rather that FITSIO.colnames so as to ignore variable length columns
Tables.getcolumn(t::EitherTableHDU, i::Int) = Tables.getcolumn(t, Tables.columnnames(t)[i])

end
