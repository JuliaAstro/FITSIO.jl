using WCS: WCSTransform
using FITSIO: fitswrite, read_header
using Printf: @sprintf

norm(text) = replace(text, r"'(-?\d+\.?\d*)\s*'" => m -> begin
    num = parse(Float64, match(r"-?\d+\.?\d*", m).match)
    @sprintf "%.3f" num
end)

@testset "WCS handling" begin
    # Create sample fits data
    img = [6 7; 8 9]
    wcs = WCSTransform(2;
        cdelt = [-0.066667, 0.066667],
        ctype = ["RA---AIR", "DEC--AIR"],
        crpix = [-234.75, 8.3393],
        crval = [0., -90],
        pv    = [(2, 1, 45.0)],
    )
    header_wcs = FITSHeader(wcs)

    # Check output
    header_default_str = """SIMPLE  =                    T / file does conform to FITS standard
       BITPIX  =                   64 / number of bits per data pixel
       NAXIS   =                    2 / number of data axes
       NAXIS1  =                    2 / length of data axis 1
       NAXIS2  =                    2 / length of data axis 2
       EXTEND  =                    T / FITS dataset may contain extensions
       COMMENT   FITS (Flexible Image Transport System) format is defined in 'Astronom
       COMMENT   and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H
       """

    header_wcs_str = """WCSAXES = '2       '           / Number of coordinate axes
       CRPIX1  = '-234.7500'          / Pixel coordinate of reference point
       CRPIX2  = '8.3393  '           / Pixel coordinate of reference point
       CDELT1  = '-0.066667'          / [deg] Coordinate increment at reference point
       CDELT2  = '0.066667'           / [deg] Coordinate increment at reference point
       CUNIT1  = 'deg     '           / Units of coordinate increment and value
       CUNIT2  = 'deg     '           / Units of coordinate increment and value
       CTYPE1  = 'RA---AIR'           / Right ascension, Airys zenithal projection
       CTYPE2  = 'DEC--AIR'           / Declination, Airys zenithal projection
       CRVAL1  = '0.0     '           / [deg] Coordinate value at reference point
       CRVAL2  = '-90.0   '           / [deg] Coordinate value at reference point
       PV2_1   = '45.0    '           / AIR projection parameter
       LONPOLE = '180.0   '           / [deg] Native longitude of celestial pole
       LATPOLE = '-90.0   '           / [deg] Native latitude of celestial pole
       MJDREF  = '0.0     '           / [d] MJD of fiducial time
       RADESYS = 'ICRS    '           / Equatorial coordinate system
       COMMENT WCS header keyrecords produced by WCSLIB 7.7"""

    @test (norm ∘ string)(header_wcs) == norm(header_wcs_str)

    tempnamefits() do fname
        fitswrite(fname, img; header = header_wcs)
        @test (norm ∘ string ∘ read_header)(fname) == norm(header_default_str * header_wcs_str)
    end
end
