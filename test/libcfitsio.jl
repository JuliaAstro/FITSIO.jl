using FITSIO.Libcfitsio

function writehealpix(filename, pixels, nside, ordering, coordsys)
    if eltype(pixels) == Float32
        tform = "1E"
    elseif eltype(pixels) == Float64
        tform = "1D"
    end 

    file = fits_create_file("!"*filename)
    try
        fits_create_img(file, Int16, Int[])
        fits_write_date(file)
        fits_movabs_hdu(file, 1)
        fits_create_binary_tbl(file, length(pixels), [("SIGNAL", tform, "")], "BINTABLE")
        fits_write_key(file, "PIXTYPE",  "HEALPIX", "HEALPIX pixelization")
        fits_write_key(file, "ORDERING", ordering,  "Pixel ordering scheme (either RING or NESTED)")
        fits_write_key(file, "NSIDE",    nside,     "Resolution parameter for HEALPIX")
        fits_write_key(file, "COORDSYS", coordsys,  "Pixelization coordinate system")
        fits_write_comment(file, "G = galactic, E = ecliptic, C = celestial = equatorial")
        fits_write_col(file, 1, 1, 1, pixels)
    finally
        fits_close_file(file)
    end 
end

function readhealpix(filename)
    file = fits_open_file(filename)
    try
        hdutype = fits_movabs_hdu(file, 2)
        tform, tform_comment = fits_read_key_str(file, "TFORM1")
        if tform == "1E"
            T = Float32
        elseif tform == "1D"
            T = Float64
        end

        naxes, naxis_comment = fits_read_key_lng(file, "NAXIS")
        naxis, nfound = fits_read_keys_lng(file, "NAXIS", 1, naxes)
        nside, nside_comment = fits_read_key_lng(file, "NSIDE")
        npix = 12*nside*nside

        ordering, ordering_comment = fits_read_key_str(file, "ORDERING")
        coordsys, coordsys_comment = fits_read_key_str(file, "COORDSYS")

        pixels = zeros(T, npix)
        fits_read_col(file, 1, 1, 1, pixels)

        return pixels, nside, ordering, coordsys

    finally
        fits_close_file(file)
    end
end


# test reading/writing Healpix maps as FITS binary tables using the Libcfitsio interface
for T in (Float32, Float64)
    nside = 4
    npix = 12*nside*nside
    pixels = rand(T, npix)
    ordering = "NESTED"
    coordsys = "G"

    filename = tempname() * ".fits"
    try
        writehealpix(filename, pixels, nside, ordering, coordsys)
        @test readhealpix(filename) == (pixels, nside, ordering, coordsys)
    finally
        ispath(filename) && rm(filename)
    end
end

