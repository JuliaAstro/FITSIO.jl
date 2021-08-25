# FITS methods

const VERBOSE_MODE = Dict("r"=>"read-only",
                          "w"=>"read-write",
                          "r+"=>"append")

# helper function for show()
function show_ascii_table(io, names, cols, spaces=2, indent=0)
    ncols = length(cols)
    ncols >= 1 || error("No columns")
    nrows = length(cols[1])
    length(names) == ncols || error("length of cols and names must match")
    for i=1:ncols
        length(cols[i]) == nrows || error("column length mismatch")
    end

    lengths = [max(maximum(length, cols[i]), length(names[i])) + spaces
               for i=1:ncols]
    for i = 1:ncols
        print(io, rpad(names[i], lengths[i]))
    end
    println(io)
    for j = 1:nrows
        print(io, " "^indent)
        for i=1:ncols
            print(io, rpad(cols[i][j], lengths[i]))
        end
        j != nrows && println(io)
    end
end

"""
    length(f::FITS)

Number of HDUs in the file.
"""
function length(f::FITS)
    fits_assert_open(f.fitsfile)
    Int(fits_get_num_hdus(f.fitsfile))
end

lastindex(f::FITS) = length(f)

# Iteration
iterate(f::FITS, state=1) =
    (state ≤ length(f) ? (f[state], state + 1) : nothing)

function show(io::IO, f::FITS)
    fits_assert_open(f.fitsfile)

    print(io, """File: $(f.filename)
    Mode: $(repr(f.mode)) ($(VERBOSE_MODE[f.mode]))
    """)

    nhdu = length(f)

    if nhdu == 0
        print(io, "No HDUs.")
    else
        print(io, "HDUs: ")

        names = Vector{String}(undef, nhdu)
        vers  = Vector{String}(undef, nhdu)
        types = Vector{String}(undef, nhdu)
        for i = 1:nhdu
            t = fits_movabs_hdu(f.fitsfile, i)
            types[i] = (t == :image_hdu ? "Image" :
                        t == :binary_table ? "Table" :
                        t == :ascii_table ? "ASCIITable" :
                        error("unknown HDU type"))
            names[i] = something(fits_try_read_extname(f.fitsfile), "")
            ver = fits_try_read_extver(f.fitsfile)
            vers[i] = ver === nothing ? "" : string(ver)
        end

        nums = [string(i) for i=1:nhdu]

        # only display version info if present
        if maximum(length, vers) > 0
            dispnames = ["Num", "Name", "Ver", "Type"]
            dispcols = Vector{String}[nums, names, vers, types]
        else
            dispnames = ["Num", "Name", "Type"]
            dispcols = Vector{String}[nums, names, types]
        end

        show_ascii_table(io, dispnames, dispcols, 2, 6)
    end
end

# Returns HDU object based on extension number
function getindex(f::FITS, i::Integer)
    fits_assert_open(f.fitsfile)

    if haskey(f.hdus, i)
        return f.hdus[i]
    end

    if i > length(f)
        error("index out of bounds")
    end
    hdutype = fits_movabs_hdu(f.fitsfile, i)
    f.hdus[i] = (hdutype == :image_hdu ? ImageHDU(f.fitsfile, i) :
                 hdutype == :binary_table ? TableHDU(f.fitsfile, i) :
                 hdutype == :ascii_table ? ASCIITableHDU(f.fitsfile, i) :
                 error("bad HDU type"))
    return f.hdus[i]
end

# Returns HDU based on hduname, version
function getindex(f::FITS, name::AbstractString, ver::Int=0)
    fits_assert_open(f.fitsfile)
    fits_movnam_hdu(f.fitsfile, name, ver)
    i = fits_get_hdu_num(f.fitsfile)

    if haskey(f.hdus, i)
        return f.hdus[i]
    end

    hdutype = fits_get_hdu_type(f.fitsfile)
    f.hdus[i] = (hdutype == :image_hdu ? ImageHDU(f.fitsfile, i) :
                 hdutype == :binary_table ? TableHDU(f.fitsfile, i) :
                 hdutype == :ascii_table ? ASCIITableHDU(f.fitsfile, i) :
                 error("bad HDU type"))
    return f.hdus[i]
end

"""
    deleteat!(f::FITS, i::Integer)

Delete the HDU at index `i` in the FITS file. If `i == 1`, this deletes the primary HDU and replaces it
with a bare HDU with no data and a minimal header. If `i > 1`, this removes the HDU at index `i`
and moves the following HDUs forward.
"""
function deleteat!(f::FITS, i::Integer)
    fits_assert_open(f.fitsfile)

    hdu = f[i]
    hdu.ext = -1 # indicate that the hdu is deleted
    delete!(f.hdus, i)
    fits_movabs_hdu(f.fitsfile, i)
    fits_delete_hdu(f.fitsfile)
    f
end

isdeleted(hdu) = hdu.ext == -1
assert_exists(hdu) = isdeleted(hdu) && error("HDU doesn't exist, it has been deleted previously")
assert_open(hdu) = assert_exists(hdu) && fits_assert_open(hdu.fitsfile)

"""
    close(f::FITS)

Close the file.

Subsequent attempts to operate on `f` will result in an error. `FITS` objects are also
automatically closed when they are garbage collected.
"""
function close(f::FITS)
    fits_assert_open(f.fitsfile)
    fits_close_file(f.fitsfile)
    f.filename = ""
    f.mode = ""
    empty!(f.hdus)
    nothing
end
