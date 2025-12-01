module WCSExt

using WCS: WCSTransform, to_header
import FITSIO: FITSHeader

function FITSHeader(wcs::WCSTransform)
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
