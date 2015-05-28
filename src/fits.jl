# FITS methods

function length(f::FITS)
    fits_assert_open(f.fitsfile)
    @compat Int(fits_get_num_hdus(f.fitsfile))
end

endof(f::FITS) = length(f)

function show(io::IO, f::FITS)
    fits_assert_open(f.fitsfile)

    # Get name and type of all HDUs.
    nhdu = length(f)
    names = Array(ASCIIString, nhdu)
    vers = Array(ASCIIString, nhdu)
    types = Array(ASCIIString, nhdu)
    for i = 1:nhdu
        t = fits_movabs_hdu(f.fitsfile, i)
        types[i] = (t == :image_hdu ? "Image" :
                    t == :binary_table ? "Table" :
                    t == :ascii_table ? "ASCIITable" :
                    error("unknown HDU type"))
        nname = _try_read_key(f.fitsfile, ASCIIString, ("EXTNAME", "HDUNAME"))
        names[i] = get(nname, "")
        nver = _try_read_key(f.fitsfile, Int, ("EXTVER", "HDUVER"))
        vers[i] = isnull(nver) ? "" : string(get(nver))
    end

    namelen = max(maximum(length, names), 4) + 2
    verlen = maximum(length, vers)
    if verlen > 0
        verlen = max(verlen, 3) + 2
        namehead = string(rpad("Name", namelen), rpad("Ver", verlen))
    else
        namehead = rpad("Name", namelen)
    end

    print(io, """File: $(f.filename)
    Mode: \"$(f.mode)\"
    HDUs: Num   $(namehead)Type
    """)
    for i in 1:nhdu
        @printf io "      %-5d %s%s%s\n" i rpad(names[i], namelen) rpad(vers[i], verlen) types[i]
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
function getindex(f::FITS, name::String, ver::Int=0)
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
