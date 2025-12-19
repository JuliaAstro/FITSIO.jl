module WCSExt

using WCS: WCSTransform, to_header
import FITSIO: FITSIO, FITSHeader

"""
    FITSHeader(wcs::WCS.WCSTransform)

Construct a [`FITSHeader`](@ref) from a [`WCSTransform`](@extref WCS.WCSTransform) supplied by [WCS.jl](@extref).

# Examples

```jldoctest
julia> using FITSIO, WCS

julia> wcs = WCSTransform(2;
           cdelt = [-0.066667, 0.066667],
           ctype = ["RA---AIR", "DEC--AIR"],
           crpix = [-234.75, 8.3393],
           crval = [0., -90],
           pv    = [(2, 1, 45.0)],
       )
WCSTransform(naxis=2, cdelt=[-0.066667, 0.066667], crval=[0.0, -90.0], crpix=[-234.75, 8.3393])

julia> FITSHeader(wcs)
WCSAXES = '2       '           / Number of coordinate axes
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
COMMENT WCS header keyrecords produced by WCSLIB 7.7
```
"""
function FITSIO.FITSHeader(wcs::WCSTransform)
	# Split string into 80-character card images
	card_images = Iterators.partition(to_header(wcs), 80)

	# Remove any blank lines
	is_empty = isempty ∘ strip
	card_images = Iterators.filter(card_images) do card_image
		!is_empty(card_image)
	end

	# Split each of those card images into their (key, value, comment) parts
	card_image_parts = map(card_images) do card_image
        card_image = replace(card_image, "'" => "") # Remove single quotes
		map(strip, split(card_image, ['=', '/' ]))
	end

	# Deal with the special comment case
    comment_card = (first ∘ pop!)(card_image_parts)
	comment = (strip ∘ last)(split(comment_card, "COMMENT"))
	push!(card_image_parts, ["COMMENT", "", comment])

	# Store
	k, v, c = eachcol(stack(card_image_parts; dims = 1))
	return FITSHeader(string.(k), string.(v), string.(c))
end

end # module
