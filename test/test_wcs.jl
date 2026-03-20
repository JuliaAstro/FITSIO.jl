using WCS: WCSTransform
using FITSIO: fitswrite, read_header

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

    @test header_wcs isa FITSHeader

    # Verify that numeric WCS keywords are stored with the correct types,
    # not as strings (per FITS standard compliance).
    @test header_wcs["WCSAXES"] isa Int
    @test header_wcs["WCSAXES"] == 2

    @test header_wcs["CRVAL1"] isa Float64
    @test header_wcs["CRVAL1"] == 0.0
    @test header_wcs["CRVAL2"] isa Float64
    @test header_wcs["CRVAL2"] == -90.0

    @test header_wcs["CRPIX1"] isa Float64
    @test header_wcs["CRPIX1"] == -234.75
    @test header_wcs["CRPIX2"] isa Float64
    @test header_wcs["CRPIX2"] == 8.3393

    @test header_wcs["CDELT1"] isa Float64
    @test header_wcs["CDELT1"] == -0.066667
    @test header_wcs["CDELT2"] isa Float64
    @test header_wcs["CDELT2"] == 0.066667

    # String-valued WCS keywords should remain as strings
    @test header_wcs["CTYPE1"] isa String
    @test header_wcs["CTYPE1"] == "RA---AIR"
    @test header_wcs["CTYPE2"] isa String
    @test header_wcs["CTYPE2"] == "DEC--AIR"
    @test header_wcs["RADESYS"] isa String
    @test header_wcs["RADESYS"] == "ICRS"

    # Round-trip test: write to file and read back, verify types are preserved
    tempnamefits() do fname
        fitswrite(fname, img; header = header_wcs)
        hdr_roundtrip = read_header(fname)
        @test hdr_roundtrip["WCSAXES"] isa Int
        @test hdr_roundtrip["WCSAXES"] == 2
        @test hdr_roundtrip["CRVAL1"] isa Float64
        @test hdr_roundtrip["CRVAL1"] == 0.0
        @test hdr_roundtrip["CRVAL2"] isa Float64
        @test hdr_roundtrip["CRVAL2"] == -90.0
        @test hdr_roundtrip["CTYPE1"] isa String
        @test hdr_roundtrip["CTYPE1"] == "RA---AIR"
    end
end
