# FITS methods

const VERBOSE_MODE = @compat Dict("r"=>"read-only",
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
    print(io, "\n")
    for j = 1:nrows
        print(io, " "^indent)
        for i=1:ncols
            print(io, rpad(cols[i][j], lengths[i]))
        end
        print(io, "\n")
    end
end

function length(f::FITS)
    fits_assert_open(f.fitsfile)
    @compat Int(fits_get_num_hdus(f.fitsfile))
end

endof(f::FITS) = length(f)

# Iteration
start(f::FITS) = 1
next(f::FITS, state) = (f[state], state + 1)
done(f::FITS, state) = state > length(f)

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

        names = Array(Compat.ASCIIString, nhdu)
        vers = Array(Compat.ASCIIString, nhdu)
        types = Array(Compat.ASCIIString, nhdu)
        for i = 1:nhdu
            t = fits_movabs_hdu(f.fitsfile, i)
            types[i] = (t == :image_hdu ? "Image" :
                        t == :binary_table ? "Table" :
                        t == :ascii_table ? "ASCIITable" :
                        error("unknown HDU type"))
            nname = fits_try_read_extname(f.fitsfile)
            names[i] = get(nname, "")
            nver = fits_try_read_extver(f.fitsfile)
            vers[i] = isnull(nver) ? "" : string(get(nver))
        end
        
        nums = [string(i) for i=1:nhdu]

        # only display version info if present
        if maximum(length, vers) > 0
            dispnames = ["Num", "Name", "Ver", "Type"]
            dispcols = Vector{Compat.ASCIIString}[nums, names, vers, types]
        else
            dispnames = ["Num", "Name", "Type"]
            dispcols = Vector{Compat.ASCIIString}[nums, names, types]
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


function close(f::FITS)
    fits_assert_open(f.fitsfile)
    fits_close_file(f.fitsfile)
    f.filename = ""
    f.mode = ""
    empty!(f.hdus)
    nothing
end
