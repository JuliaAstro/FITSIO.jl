var documenterSearchIndex = {"docs":
[{"location":"api/#API-Reference","page":"API Reference","title":"API Reference","text":"","category":"section"},{"location":"api/#File-operations","page":"API Reference","title":"File operations","text":"","category":"section"},{"location":"api/","page":"API Reference","title":"API Reference","text":"FITS\nlength\nclose\ndeleteat!","category":"page"},{"location":"api/#FITSIO.FITS","page":"API Reference","title":"FITSIO.FITS","text":"FITS(filename::String[, mode::String = \"r\"]; extendedparser = true)\n\nOpen or create a FITS file. mode can be one of \"r\" (read-only), \"r+\" (read-write) or \"w\" (write). In \"write\" mode, any existing file of the same name is overwritten.\n\nA FITS object is a collection of \"Header-Data Units\" (HDUs) and supports the following operations:\n\nf[i]: Return the i-th HDU.\nf[name] or f[name, ver]: Return the HDU containing the given the given EXTNAME (or HDUNAME) keyword (a String), and optionally the given EXTVER (or HDUVER) number (an Integer).\nIteration:\nfor hdu in f\n    ...\nend\n\nThe keyword argument extendedparser may be used to enable or disable the extended filename parser. If disabled, filename is treated exactly as the name of the file and is not tokenized into parameters.\n\n\n\n\n\n","category":"type"},{"location":"api/#Base.length","page":"API Reference","title":"Base.length","text":"length(f::FITS)\n\nNumber of HDUs in the file.\n\n\n\n\n\nlength(hdr::FITSHeader)\n\nNumber of records in header of HDU.\n\n\n\n\n\nlength(hdu::ImageHDU)\n\nGet total number of pixels in image (product of size(hdu)).\n\n\n\n\n\n","category":"function"},{"location":"api/#Base.close","page":"API Reference","title":"Base.close","text":"close(f::FITS)\n\nClose the file.\n\nSubsequent attempts to operate on f will result in an error. FITS objects are also automatically closed when they are garbage collected.\n\n\n\n\n\n","category":"function"},{"location":"api/#Base.deleteat!","page":"API Reference","title":"Base.deleteat!","text":"deleteat!(f::FITS, i::Integer)\n\nDelete the HDU at index i in the FITS file. If i == 1, this deletes the primary HDU and replaces it with a bare HDU with no data and a minimal header. If i > 1, this removes the HDU at index i and moves the following HDUs forward.\n\n\n\n\n\n","category":"function"},{"location":"api/#Header-operations","page":"API Reference","title":"Header operations","text":"","category":"section"},{"location":"api/","page":"API Reference","title":"API Reference","text":"read_key\nwrite_key\nread_header\nFITSHeader\nlength(::FITSHeader)\nhaskey(::FITSHeader, ::String)\nkeys(::FITSHeader)\nvalues(::FITSHeader)\nget_comment\nset_comment!\ndefault_header","category":"page"},{"location":"api/#FITSIO.read_key","page":"API Reference","title":"FITSIO.read_key","text":"read_key(hdu::HDU, key::String) -> (value, comment)\n\nRead the HDU header record specified by keyword and return a tuple where value is the keyword parsed value (of type String, Bool, Int, Float64 or Nothing), comment is the keyword comment (as a string). Throw an error if key is not found.\n\nread_key(hdu::HDU, key::Integer) -> (keyname, value, comment)\n\nSame as above but FITS card is specified by its position and returns a 3 element tuple where keyname is the keyword name (a string).\n\n\n\n\n\n","category":"function"},{"location":"api/#FITSIO.write_key","page":"API Reference","title":"FITSIO.write_key","text":"write_key(hdu::HDU, key::String, value[, comment])\n\nWrite a keyword value the HDU's header. value can be a standard header type (String, Bool, Integer, AbstractFloat) or nothing, in which case the value part of the record will be empty. If the keyword already exists, the value will be overwritten. The comment will only be overwritten if given. If the keyword does not already exist, a new record will be appended at the end of the header.\n\n\n\n\n\n","category":"function"},{"location":"api/#FITSIO.read_header","page":"API Reference","title":"FITSIO.read_header","text":"read_header(filename::AbstractString, hduindex = 1) -> FITSHeader\n\nConvenience function to read the entire header corresponding to the HDU at index hduindex contained in the FITS file named filename. Functionally read_header(filename, hduindex) is equivalent to\n\nFITS(filename, \"r\") do f\n    read_header(f[hduindex])\nend\n\n\n\n\n\nread_header(hdu::HDU) -> FITSHeader\n\nRead the entire header from the given HDU and return a FITSHeader object. The value of each header record is parsed as Int, Float64, String, Bool or nothing according to the FITS standard.\n\nIf the value cannot be parsed according to the FITS standard, the value is stored as the raw unparsed String.\n\n\n\n\n\nread_header(hdu::HDU, String) -> String\n\nRead the entire header from the given HDU as a single string.\n\n\n\n\n\n","category":"function"},{"location":"api/#FITSIO.FITSHeader","page":"API Reference","title":"FITSIO.FITSHeader","text":"FITSHeader(keys::Vector{String}, values::Vector, comments::Vector{String})\n\nAn in-memory representation of the header of an HDU. It stores the (key, value, comment) information for each 80-character \"card\" in a header.\n\nNote that this structure is not linked to a FITS file in any way; it is just a convenient structure for storing the header contents after reading from a file. (This is similar to how an Array returned by read(f[1]) is not linked to the FITS file f.)  Manipulating a FITSHeader will therefore have no immediate impact on any file, even if it was created by read_header(::HDU).  You can, however, write a FITSHeader to a file using the write(::FITS, ...) methods that append a new HDU to a file.\n\n\n\n\n\n","category":"type"},{"location":"api/#Base.length-Tuple{FITSHeader}","page":"API Reference","title":"Base.length","text":"length(hdr::FITSHeader)\n\nNumber of records in header of HDU.\n\n\n\n\n\n","category":"method"},{"location":"api/#Base.haskey-Tuple{FITSHeader, String}","page":"API Reference","title":"Base.haskey","text":"haskey(hdr::FITSHeader, key::String)\n\nReturns true if key exists in header, otherwise false.\n\n\n\n\n\n","category":"method"},{"location":"api/#Base.keys-Tuple{FITSHeader}","page":"API Reference","title":"Base.keys","text":"keys(hdr::FITSHeader)\n\nArray of keywords in header of HDU (not a copy).\n\n\n\n\n\n","category":"method"},{"location":"api/#Base.values-Tuple{FITSHeader}","page":"API Reference","title":"Base.values","text":"values(hdr::FITSHeader)\n\nArray of values in header of HDU (not a copy).\n\n\n\n\n\n","category":"method"},{"location":"api/#FITSIO.get_comment","page":"API Reference","title":"FITSIO.get_comment","text":"get_comment(hdr::FITSHeader, key_or_index::Union{String,Integer})\n\nGet the comment based on keyword or index.\n\n\n\n\n\n","category":"function"},{"location":"api/#FITSIO.set_comment!","page":"API Reference","title":"FITSIO.set_comment!","text":"set_comment!(hdr::FITSHeader, key_or_index::Union{String,Integer}, comment::String)\n\nSet the comment based on keyword or index.\n\n\n\n\n\n","category":"function"},{"location":"api/#FITSIO.default_header","page":"API Reference","title":"FITSIO.default_header","text":"default_header(data::AbstractArray)\n\nCreates a default header for the given array with the SIMPLE, BITPIX, NAXIS, NAXIS*, and EXTEND entries.\n\n\n\n\n\n","category":"function"},{"location":"api/#Image-operations","page":"API Reference","title":"Image operations","text":"","category":"section"},{"location":"api/","page":"API Reference","title":"API Reference","text":"read(::ImageHDU)\nread!\nFITSIO.fitsread\nwrite(::FITS, ::StridedArray{<:Real})\nwrite(::ImageHDU, ::StridedArray{<:Real})\nFITSIO.fitswrite\neltype(::ImageHDU)\nndims(::ImageHDU)\nsize(::ImageHDU)\nlength(::ImageHDU)\ncopy_section","category":"page"},{"location":"api/#Base.read-Tuple{ImageHDU}","page":"API Reference","title":"Base.read","text":"read(hdu::ImageHDU)\nread(hdu::ImageHDU, range...)\n\nRead the data array or a subset thereof from disk. The first form reads the entire data array. The second form reads a slice of the array given by the specified ranges or integers. Dimensions specified by integers will be dropped in the returned array, while those specified by ranges will be retained.\n\nnote: Note\nJulia follows a column-major array indexing convention, so the indices provided must account for this. In particular this means that FITS files created externally following a row-major convention (eg. using astropy) will have the sequence of axes flipped when read in using FITSIO.\n\n\n\n\n\n","category":"method"},{"location":"api/#Base.read!","page":"API Reference","title":"Base.read!","text":"read!(hdu::ImageHDU, A::StridedArray)\nread!(hdu::ImageHDU, A::StridedArray, range...)\n\nRead the data or a subset thereof from disk, and save it in a pre-allocated output array A. The first form reads the entire data from disk. The second form reads a slice of the array given by the specified ranges or integers. The array A needs to have the same length as the number of elements to be read in. Additionally A needs to be stored contiguously in memory.\n\nnote: Note\nJulia follows a column-major array indexing convention, so the indices provided must account for this. In particular this means that FITS files created externally following a row-major convention (eg. using astropy) will have the sequence of the axes flipped when read in using FITSIO.\n\n\n\n\n\n","category":"function"},{"location":"api/#FITSIO.fitsread","page":"API Reference","title":"FITSIO.fitsread","text":"fitsread(filename::AbstractString[, hduindex = 1[, arrayindices...]]; extendedparser = true)\n\nConvenience function to read in an image corresponding to the HDU at index hduindex contained in the FITS file named filename. If arrayindices are provided, only a slice of the image corresponding to the indices is read in.\n\nFunctionally fitsread(filename, hduindex, arrayindices...; extendedparser) is equivalent to\n\nFITS(filename, \"r\"; extendedparser = extendedparser) do f\n    read(f[hduindex], arrayindices...)\nend\n\nThe keyword argument extendedparser may be used to enable or disable the extended filename parser. If disabled, filename is treated exactly as the name of the file and is not tokenized into parameters.\n\nnote: Note\nJulia follows a column-major array indexing convention, so the indices provided must account for this. In particular this means that FITS files created externally following a row-major convention (eg. using astropy) will have the sequence of axes flipped when read in using FITSIO.\n\nSee also: read\n\n\n\n\n\n","category":"function"},{"location":"api/#Base.write-Tuple{FITS, StridedArray{<:Real}}","page":"API Reference","title":"Base.write","text":"write(f::FITS, data::StridedArray{<:Real}; header=nothing, name=nothing, ver=nothing)\n\nAdd a new image HDU to FITS file f with contents data. The following array element types are supported: UInt8, Int8, UInt16, Int16, UInt32, Int32, Int64, Float32, Float64. If a FITSHeader object is passed as the header keyword argument, the header will also be added to the new HDU. The data to be written out must be stored contiguously in memory.\n\ntip: Unsupported element types\nIt might be possible to write out an array with an element type other than those mentioned above by reinterpreting it as one that is supported. For example, to write out a Complex array and read it back in, we may usejulia> a = rand(ComplexF64, 2)\n2-element Array{Complex{Float64},1}:\n 0.4943325325752195 + 0.2034650017475852im\n 0.2495752009567498 + 0.819163869249041im\n\n# We may write this out as Float64\njulia> FITSIO.fitswrite(\"temp.fits\", reinterpret(Float64, a))\n\n# reinterpret it back as a complex one while reading it in\njulia> reinterpret(ComplexF64, FITSIO.fitsread(\"temp.fits\"))\n2-element reinterpret(Complex{Float64}, ::Array{Float64,1}):\n 0.4943325325752195 + 0.2034650017475852im\n 0.2495752009567498 + 0.819163869249041imWhile this often works in practice, such a workaround is not officially supported by FITSIO, and care must be taken to ensure the correctness of data.\n\n\n\n\n\n","category":"method"},{"location":"api/#Base.write-Tuple{ImageHDU, StridedArray{<:Real}}","page":"API Reference","title":"Base.write","text":"write(hdu::ImageHDU, data::StridedArray{<:Real})\n\nWrite data to an existing image HDU. The data to be written out must be stored contiguously in memory.\n\n\n\n\n\n","category":"method"},{"location":"api/#FITSIO.fitswrite","page":"API Reference","title":"FITSIO.fitswrite","text":"fitswrite(filename::AbstractString, data; extendedparser = true, kwargs...)\n\nConvenience function to write the image array data to a file named filename.\n\nFunctionally fitswrite(filename, data; extendedparser, kwargs...) is equivalent to\n\nFITS(filename, \"w\"; extendedparser = extendedparser) do f\n    write(f, data; kwargs...)\nend\n\nThe keyword argument extendedparser may be used to enable or disable the extended filename parser. If disabled, filename is treated exactly as the name of the file and is not tokenized into parameters.\n\nwarn: Warning\nExisting files with the same name will be overwritten.\n\nSee also: write\n\n\n\n\n\n","category":"function"},{"location":"api/#Base.eltype-Tuple{ImageHDU}","page":"API Reference","title":"Base.eltype","text":"eltype(hdu::ImageHDU)\n\nReturn the element type of the image in hdu.\n\n\n\n\n\n","category":"method"},{"location":"api/#Base.ndims-Tuple{ImageHDU}","page":"API Reference","title":"Base.ndims","text":"ndims(hdu::ImageHDU)\n\nGet number of image dimensions, without reading the image into memory.\n\n\n\n\n\n","category":"method"},{"location":"api/#Base.size-Tuple{ImageHDU}","page":"API Reference","title":"Base.size","text":"size(hdu::ImageHDU)\nsize(hdu::ImageHDU, i)\n\nGet image dimensions (or ith dimension), without reading the image into memory.\n\n\n\n\n\n","category":"method"},{"location":"api/#Base.length-Tuple{ImageHDU}","page":"API Reference","title":"Base.length","text":"length(hdu::ImageHDU)\n\nGet total number of pixels in image (product of size(hdu)).\n\n\n\n\n\n","category":"method"},{"location":"api/#FITSIO.copy_section","page":"API Reference","title":"FITSIO.copy_section","text":"copy_section(hdu, dest, r...)\n\nCopy a rectangular section of an image and write it to a new FITS primary image or image extension in FITS object dest. The new image HDU is appended to the end of dest. All the keywords in the input image will be copied to the output image. The common WCS keywords will be updated if necessary to correspond to the coordinates of the section.\n\nExamples\n\nCopy the lower-left 200 x 200 pixel section of the image in hdu to an open file, f\n\ncopy_section(hdu, f, 1:200, 1:200)\n\nSame as above but only copy odd columns in y:\n\ncopy_section(hdu, f, 1:200, 1:2:200)\n\n\n\n\n\n","category":"function"},{"location":"api/#Table-operations","page":"API Reference","title":"Table operations","text":"","category":"section"},{"location":"api/","page":"API Reference","title":"API Reference","text":"FITSIO.colnames\nwrite(::FITS, ::Dict{String})\nwrite(::FITS, ::Vector{String}, ::Vector)\nread(::TableHDU, ::String)","category":"page"},{"location":"api/#FITSIO.colnames","page":"API Reference","title":"FITSIO.colnames","text":"colnames(hdu) -> Vector{String}\n\nReturn the names of columns in a table HDU.\n\n\n\n\n\n","category":"function"},{"location":"api/#Base.write-Tuple{FITS, Dict{String}}","page":"API Reference","title":"Base.write","text":"write(f::FITS, data::Dict; hdutype=TableHDU, name=nothing, ver=nothing, header=nothing, units=nothing, varcols=nothing)\n\nCreate a new table extension and write data to it. If the FITS file is currently empty then a dummy primary array will be created before appending the table extension to it. data should be a dictionary with String keys (giving the column names) and Array values (giving data to write to each column). The following types are supported in binary tables: UInt8, Int8, UInt16, Int16, UInt32, Int32, Int64, Float32, Float64, Complex{Float32}, Complex{Float64}, String, Bool.\n\nOptional inputs:\n\nhdutype: Type of table extension to create. Can be either TableHDU (binary table) or ASCIITableHDU (ASCII table).\nname: Name of extension.\nver: Version of extension (Int).\nheader: FITSHeader instance to write to new extension.\nunits: Dictionary mapping column name to units (as a string).\nvarcols: An array giving the column names or column indicies to write as \"variable-length columns\".\n\nnote: Variable length columns\nVariable length columns allow a column's row entries to contain arrays of different lengths. They can potentially save diskspace when the rows of a column vary greatly in length, as the column data is all written to a contiguous heap area at the end of the table. Only column data of type Vector{String} or types such as Vector{Vector{UInt8}} can be written as variable length columns. In the second case, ensure that the column data type is a leaf type. That is, the type cannot be Vector{Vector{T}}, which would be an array of arrays having potentially non-uniform element types (which would not be writable as a FITS table column).\n\n\n\n\n\n","category":"method"},{"location":"api/#Base.write-Tuple{FITS, Vector{String}, Vector}","page":"API Reference","title":"Base.write","text":"write(f::FITS, colnames, coldata; hdutype=TableHDU, name=nothing, ver=nothing, header=nothing, units=nothing, varcols=nothing)\n\nSame as write(f::FITS, data::Dict; ...) but providing column names and column data as a separate arrays. This is useful for specifying the order of the columns. Column names must be Vector{String} and column data must be a vector of arrays.\n\n\n\n\n\n","category":"method"},{"location":"api/#Base.read-Tuple{TableHDU, String}","page":"API Reference","title":"Base.read","text":"read(hdu::TableHDU, colname; case_sensitive=true)\n\nRead a column as an array from the given table HDU.\n\nThe column name may contain wild card characters (*, ?, or #). The * wild card character matches any sequence of characters (including zero characters) and the ? character matches any single character. The # wildcard will match any consecutive string of decimal digits (0-9). The string must match a unique column.  The optional boolean keyword case_sensitive, true by default, specifies whether the column name is to be considered case sensitive.\n\nnote: Array order\nJulia arrays are column-major (like Fortran), not row-major (like C and numpy), so elements of multi-dimensional columns will be the transpose of what you get with astropy.\n\n\n\n\n\n","category":"method"},{"location":"#FITSIO.jl","page":"Introduction","title":"FITSIO.jl","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"(Image: GitHub) (Image: Build Status) (Image: Coverage Status)","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"A Julia package for reading and writing Flexible Image Transport System (FITS) files, based on the cfitsio library.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"The interface is inspired by Erin Sheldon's fitsio Python package.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"warning: Warning\nThe Libcfitsio submodule has been moved to CFITSIO.jl and will be deprecated in a future release.","category":"page"},{"location":"#Installation","page":"Introduction","title":"Installation","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"FITSIO.jl can be installed using the built-in package manager","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"pkg> add FITSIO","category":"page"},{"location":"#Usage","page":"Introduction","title":"Usage","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"To open an existing file for reading:","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"julia> using FITSIO\n\njulia> f = FITS(\"file.fits\")\nFile: file.fits\nMode: \"w\" (read-write)\nHDUs: Num  Name  Type   \n      1          Image  \n      2          Table  ","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"(At the REPL, information about the file contents is shown.)","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"A FITS file consists of one or more header-data units (HDUs), concatenated one after the other. The FITS object therefore is represented as a collection of these HDUs.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Get information about the first HDU:","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"julia> f[1]\nFile: file.fits\nHDU: 1\nType: Image\nDatatype: Float64\nDatasize: (800, 800)","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Iterate over HDUs in the file:","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"julia> for hdu in f; println(typeof(hdu)); end\nFITSIO.ImageHDU\nFITSIO.TableHDU","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Each HDU can contain image data, or table data (either binary or ASCII-formatted). For image extensions, get the size of the image without reading it:","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"julia> ndims(f[1])\n    2\n\njulia> size(f[1])\n(800,800)\n\njulia> size(f[1], 2)\n800","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Read an image from disk:","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"julia> data = read(f[1]);  # read an image from disk\n\njulia> data = read(f[1], :, 790:end);  # read just a subset of image","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Show info about a binary table:","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"julia> f[2]\nFile: file.fits\nHDU: 2\nType: Table\nRows: 20\nColumns: Name  Size  Type    TFORM  \n         col2        String  5A     \n         col1        Int64   1K     ","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Read a column from the table:","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":" julia> data = read(f[2], \"col1\")","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Table HDUs implement the Tables.jl interface, so you can load them into other table types, like DataFrames.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"julia> df = DataFrame(f[2])","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Variable length columns are not supported by the Tables.jl interface, and Tables methods will ignore them.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Read the entire header into memory and get values from it:","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"julia> header = read_header(f[1]);  # read the entire header from disk\n\njulia> length(header)  # total number of records in header\n17\n\njulia> haskey(header, \"NAXIS1\")  # check if a key exists\ntrue\n\njulia> header[\"NAXIS1\"]  # get value by keyword\n800\n\njulia> header[4]  # get value by position\n800\n\njulia> get_comment(header, \"NAXIS\")  # get comment for a given keyword\n\"length of data axis 1\"","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Read just a single header record without reading the entire header:","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"julia> read_key(f[1], 4)  # by position\n(\"NAXIS1\",800,\"length of data axis 1\")\n\njulia> read_key(f[1], \"NAXIS1\")  # read by keyword\n(800,\"length of data axis 1\")","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Manipulate a header in memory:","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"julia> header[\"NEWKEY\"] = 10  # change or add a keyword\n\njulia> set_comment!(header, \"NEWKEY\", \"this is a comment\")","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Close the file:","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"julia> close(f)","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"(FITS objects are also closed automatically when garbage collected.)","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Open a new file for writing:","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"julia> f = FITS(\"newfile.fits\", \"w\");","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"The second argument can be \"r\" (read-only; default), \"r+\" (read-write) or \"w\" (write). In \"write\" mode, any existing file of the same name is overwritten.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Write an image to the file:","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"julia> data = reshape([1:100;], 5, 20)\n\njulia> write(f, data)  # Write a new image extension with the data\njulia> close(f)","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"To write some header keywords in the new extension, pass a FITSHeader instance as a keyword: write(f, data; header=header)","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Overwrite image data in an existing file:","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"julia> f = FITS(\"newfile.fits\", \"r+\")  # Reopen the file in read-write mode\njulia> data = reshape([101:200;], 5, 20)  # Prepare new image data\njulia> image_hdu = f[1]\njulia> write(image_hdu, data)  # Overwrite the image","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Write a table to the file:","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"julia> data = Dict(\"col1\"=>[1., 2., 3.], \"col2\"=>[1, 2, 3]);\n\njulia> write(f, data)  # write a new binary table to a new extension","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"tip: Compressed storage\nSetting the file extension to .gz will automatically use GZIP compression and save on storage space.julia> FITS(\"abc.fits\", \"w\") do f # save the image uncompressed\n           write(f, ones(200,200))\n       end\n\njulia> filesize(\"abc.fits\")\n325440\n\njulia> FITS(\"abc.fits.gz\", \"w\") do f # save the image compressed\n            write(f, ones(200,200))\n       end\n\njulia> filesize(\"abc.fits.gz\")\n2117Alternately the compression algorithm might be specified in square brackets after the filename. Check the CFITSIO website for the details of this usage.julia> FITS(\"abc.fits[compress R 100,100]\", \"w\") do f # Rice algorithm with a 100 x 100 pixel tile size\n           write(f, ones(200,200))\n       end\n\njulia> filesize(\"abc.fits\")\n8640warn: Warn\nCompression is \"loss-less\" for images with integer pixel values, and might be lossy for floating-point images. ","category":"page"}]
}