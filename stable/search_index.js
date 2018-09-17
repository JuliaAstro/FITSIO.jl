var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#FITSIO.jl-1",
    "page": "Introduction",
    "title": "FITSIO.jl",
    "category": "section",
    "text": "A Julia package for reading and writing Flexible Image Transport System (FITS) files, based on the cfitsio library.The interface is inspired by Erin Sheldon\'s fitsio Python package."
},

{
    "location": "index.html#Installation-1",
    "page": "Introduction",
    "title": "Installation",
    "category": "section",
    "text": "FITSIO is available for Julia 0.6 and later versions, and can be installed with Julia\'s built-in package manager. In a Julia session run the commandjulia> Pkg.update()\njulia> Pkg.add(\"FITSIO\")On Linux or OS X, if it isn\'t already installed on your system, the cfitsio library is automatically downloaded and compiled (in your Julia packages directory). On Windows, a compiled dll will be downloaded."
},

{
    "location": "index.html#Usage-1",
    "page": "Introduction",
    "title": "Usage",
    "category": "section",
    "text": "To open an existing file for reading:julia> using FITSIO\n\njulia> f = FITS(\"file.fits\")\nFile: file.fits\nMode: \"w\" (read-write)\nHDUs: Num  Name  Type   \n      1          Image  \n      2          Table  (At the REPL, information about the file contents is shown.)A FITS file consists of one or more header-data units (HDUs), concatenated one after the other. The FITS object therefore is represented as a collection of these HDUs.Get information about the first HDU:julia> f[1]\nFile: file.fits\nHDU: 1\nType: Image\nDatatype: Float64\nDatasize: (800, 800)Iterate over HDUs in the file:julia> for hdu in f; println(typeof(hdu)); end\nFITSIO.ImageHDU\nFITSIO.TableHDUEach HDU can contain image data, or table data (either binary or ASCII-formatted). For image extensions, get the size of the image without reading it:julia> ndims(f[1])\n    2\n\njulia> size(f[1])\n(800,800)\n\njulia> size(f[1], 2)\n800Read an image from disk:julia> data = read(f[1]);  # read an image from disk\n\njulia> data = read(f[1], :, 790:end);  # read just a subset of imageShow info about a binary table:julia> f[2]\nFile: file.fits\nHDU: 2\nType: Table\nRows: 20\nColumns: Name  Size  Type    TFORM  \n         col2        String  5A     \n         col1        Int64   1K     Read a column from the table: julia> data = read(f[2], \"col1\")Read the entire header into memory and get values from it:julia> header = read_header(f[1]);  # read the entire header from disk\n\njulia> length(header)  # total number of records in header\n17\n\njulia> haskey(header, \"NAXIS1\")  # check if a key exists\ntrue\n\njulia> header[\"NAXIS1\"]  # get value by keyword\n800\n\njulia> header[4]  # get value by position\n800\n\njulia> get_comment(header, \"NAXIS\")  # get comment for a given keyword\n\"length of data axis 1\"Read just a single header record without reading the entire header:julia> read_key(f[1], 4)  # by position\n(\"NAXIS1\",800,\"length of data axis 1\")\n\njulia> read_key(f[1], \"NAXIS1\")  # read by keyword\n(800,\"length of data axis 1\")Manipulate a header in memory:julia> header[\"NEWKEY\"] = 10  # change or add a keyword\n\njulia> set_comment!(header, \"NEWKEY\", \"this is a comment\")Close the file:julia> close(f)(FITS objects are also closed automatically when garbage collected.)Open a new file for writing:julia> f = FITS(\"newfile.fits\", \"w\");The second argument can be \"r\" (read-only; default), \"r+\" (read-write) or \"w\" (write). In \"write\" mode, any existing file of the same name is overwritten.Write an image to the file:julia> data = reshape([1:100], 5, 20);\n\njulia> write(f, data)  # Write a new image extension with the dataTo write some header keywords in the new extension, pass a FITSHeader instance as a keyword: write(f, data; header=header)Write a table to the file:julia> data = Dict(\"col1\"=>[1., 2., 3.], \"col2\"=>[1, 2, 3]);\n\njulia> write(f, data)  # write a new binary table to a new extension"
},

{
    "location": "api.html#",
    "page": "API Reference",
    "title": "API Reference",
    "category": "page",
    "text": ""
},

{
    "location": "api.html#API-Reference-1",
    "page": "API Reference",
    "title": "API Reference",
    "category": "section",
    "text": ""
},

{
    "location": "api.html#FITSIO.FITS",
    "page": "API Reference",
    "title": "FITSIO.FITS",
    "category": "type",
    "text": "FITS(filename::String, mode::String=\"r\")\n\nOpen or create a FITS file. mode can be one of \"r\" (read-only), \"r+\" (read-write) or \"w\" (write). In \"write\" mode, any existing file of the same name is overwritten.\n\nA FITS object is a collection of \"Header-Data Units\" (HDUs) and supports the following operations:\n\nf[i]: Return the i-th HDU.\nf[name] or f[name, ver]: Return the HDU containing the given the given EXTNAME (or HDUNAME) keyword (a String), and optionally the given EXTVER (or HDUVER) number (an Integer).\nIteration:\nfor hdu in f\n    ...\nend\n\n\n\n\n\n"
},

{
    "location": "api.html#Base.length-Tuple{FITS}",
    "page": "API Reference",
    "title": "Base.length",
    "category": "method",
    "text": "length(f::FITS)\n\nNumber of HDUs in the file.\n\n\n\n\n\n"
},

{
    "location": "api.html#Base.close-Tuple{FITS}",
    "page": "API Reference",
    "title": "Base.close",
    "category": "method",
    "text": "close(f::FITS)\n\nClose the file.\n\nSubsequent attempts to operate on f will result in an error. FITS objects are also automatically closed when they are garbage collected.\n\n\n\n\n\n"
},

{
    "location": "api.html#File-operations-1",
    "page": "API Reference",
    "title": "File operations",
    "category": "section",
    "text": "FITS\nlength(::FITS)\nclose(::FITS)"
},

{
    "location": "api.html#FITSIO.read_key",
    "page": "API Reference",
    "title": "FITSIO.read_key",
    "category": "function",
    "text": "read_key(hdu, key::String) -> (value, comment)\n\nreads the HDU header record specified by keyword and returns a tuple where value is the keyword parsed value (of type String, Bool, Int, Float64 or Nothing), comment is the keyword comment (as a string). An error is thrown if key is not found.\n\nread_key(hdu, key::Integer) -> (keyname, value, comment)\n\nsame as above but FITS card is specified by its position and returns a 3 element tuple where keyname is the keyword name (a string).\n\n\n\n\n\n"
},

{
    "location": "api.html#FITSIO.write_key",
    "page": "API Reference",
    "title": "FITSIO.write_key",
    "category": "function",
    "text": "write_key(hdu, key::String, value[, comment])\n\nWrite a keyword value the HDU\'s header. value can be a standard header type (String, Bool, Integer, AbstractFloat) or nothing, in which case the value part of the record will be empty. If the keyword already exists, the value will be overwritten. The comment will only be overwritten if given. If the keyword does not already exist, a new record will be appended at the end of the header.\n\n\n\n\n\n"
},

{
    "location": "api.html#FITSIO.read_header",
    "page": "API Reference",
    "title": "FITSIO.read_header",
    "category": "function",
    "text": "read_header(hdu) -> FITSHeader\n\nRead the entire header from the given HDU and return a FITSHeader object. The value of each header record is parsed as Int, Float64, String, Bool or nothing according to the FITS standard.\n\nIf the value cannot be parsed according to the FITS standard, the value is stored as the raw unparsed String.\n\n\n\n\n\nread_header(hdu, String) -> String\n\nRead the entire header from the given HDU as a single string.\n\n\n\n\n\n"
},

{
    "location": "api.html#FITSIO.FITSHeader",
    "page": "API Reference",
    "title": "FITSIO.FITSHeader",
    "category": "type",
    "text": "FITSHeader(keys::Vector{String}, values::Vector, comments::Vector{String})\n\nAn in-memory representation of the header of an HDU. It stores the (key, value, comment) information for each 80-character \"card\" in a header.\n\nNote that this structure is not linked to a FITS file in any way; it is just a convenient structure for storing the header contents after reading from a file. (This is similar to how an Array returned by read(f[1]) is not linked to the FITS file f.)  Manipulating a FITSHeader will therefore have no immediate impact on any file, even if it was created by read_header(::HDU).  You can, however, write a FITSHeader to a file using the write(::FITS, ...) methods that append a new HDU to a file.\n\n\n\n\n\n"
},

{
    "location": "api.html#Base.length-Tuple{FITSHeader}",
    "page": "API Reference",
    "title": "Base.length",
    "category": "method",
    "text": "length(hdr)\n\nNumber of records.\n\n\n\n\n\n"
},

{
    "location": "api.html#Base.haskey-Tuple{FITSHeader,String}",
    "page": "API Reference",
    "title": "Base.haskey",
    "category": "method",
    "text": "haskey(hdr)\n\nHeader keyword exists.\n\n\n\n\n\n"
},

{
    "location": "api.html#Base.keys-Tuple{FITSHeader}",
    "page": "API Reference",
    "title": "Base.keys",
    "category": "method",
    "text": "keys(hdr)\n\nArray of keywords (not a copy).\n\n\n\n\n\n"
},

{
    "location": "api.html#Base.values-Tuple{FITSHeader}",
    "page": "API Reference",
    "title": "Base.values",
    "category": "method",
    "text": "values(hdr)\n\nArray of values (not a copy).\n\n\n\n\n\n"
},

{
    "location": "api.html#FITSIO.get_comment",
    "page": "API Reference",
    "title": "FITSIO.get_comment",
    "category": "function",
    "text": "get_comment(hdr, key)\n\nGet the comment based on keyword or index.\n\n\n\n\n\n"
},

{
    "location": "api.html#FITSIO.set_comment!",
    "page": "API Reference",
    "title": "FITSIO.set_comment!",
    "category": "function",
    "text": "set_comment!(hdr, key, comment)\n\nSet the comment baed on keyword or index.\n\n\n\n\n\n"
},

{
    "location": "api.html#Header-operations-1",
    "page": "API Reference",
    "title": "Header operations",
    "category": "section",
    "text": "read_key\nwrite_key\nread_header\nFITSHeader\nlength(::FITSHeader)\nhaskey(::FITSHeader, ::String)\nkeys(::FITSHeader)\nvalues(::FITSHeader)\nget_comment\nset_comment!"
},

{
    "location": "api.html#Base.write-Union{Tuple{T}, Tuple{FITS,Array{T,N} where N}} where T",
    "page": "API Reference",
    "title": "Base.write",
    "category": "method",
    "text": "write(f::FITS, data::Array; header=nothing, name=nothing, ver=nothing)\n\nAdd a new image HDU to FITS file f with contents data. The following array element types are supported: UInt8, Int8, UInt16, Int16, UInt32, Int32, Int64, Float32, Float64. If a FITSHeader object is passed as the header keyword argument, the header will also be added to the new HDU.\n\n\n\n\n\n"
},

{
    "location": "api.html#Base.read-Tuple{ImageHDU}",
    "page": "API Reference",
    "title": "Base.read",
    "category": "method",
    "text": "read(hdu::ImageHDU)\nread(hdu::ImageHDU, range...)\n\nRead the data array or a subset thereof from disk. The first form reads the entire data array. The second form reads a slice of the array given by the specified ranges or integers.\n\n\n\n\n\n"
},

{
    "location": "api.html#Base.ndims-Tuple{ImageHDU}",
    "page": "API Reference",
    "title": "Base.ndims",
    "category": "method",
    "text": "ndims(hdu::ImageHDU)\n\nGet number of image dimensions, without reading the image into memory.\n\n\n\n\n\n"
},

{
    "location": "api.html#Base.size-Tuple{ImageHDU}",
    "page": "API Reference",
    "title": "Base.size",
    "category": "method",
    "text": "size(hdu::ImageHDU)\nsize(hdu::ImageHDU, i)\n\nGet image dimensions (or ith dimension), without reading the image into memory.\n\n\n\n\n\n"
},

{
    "location": "api.html#Base.length-Tuple{ImageHDU}",
    "page": "API Reference",
    "title": "Base.length",
    "category": "method",
    "text": "length(hdu::ImageHDU)\n\nGet total number of pixels in image (product of size(hdu)).\n\n\n\n\n\n"
},

{
    "location": "api.html#FITSIO.copy_section",
    "page": "API Reference",
    "title": "FITSIO.copy_section",
    "category": "function",
    "text": "copy_section(hdu, dest, r...)\n\nCopy a rectangular section of an image and write it to a new FITS primary image or image extension in FITS object dest. The new image HDU is appended to the end of dest. All the keywords in the input image will be copied to the output image. The common WCS keywords will be updated if necessary to correspond to the coordinates of the section.\n\nExamples\n\nCopy the lower-left 200 x 200 pixel section of the image in hdu to an open file, f\n\ncopy_section(hdu, f, 1:200, 1:200)\n\nSame as above but only copy odd columns in y:\n\ncopy_section(hdu, f, 1:200, 1:2:200)\n\n\n\n\n\n"
},

{
    "location": "api.html#Image-operations-1",
    "page": "API Reference",
    "title": "Image operations",
    "category": "section",
    "text": "write{T}(::FITS, ::Array{T})\nread(::ImageHDU)\nndims(::ImageHDU)\nsize(::ImageHDU)\nlength(::ImageHDU)\ncopy_section"
},

{
    "location": "api.html#FITSIO.colnames",
    "page": "API Reference",
    "title": "FITSIO.colnames",
    "category": "function",
    "text": "colnames(hdu) -> Vector{String}\n\nReturn the names of columns in a table HDU.\n\n\n\n\n\n"
},

{
    "location": "api.html#Base.write-Tuple{FITS,Dict{String,V} where V}",
    "page": "API Reference",
    "title": "Base.write",
    "category": "method",
    "text": "write(f::FITS, data::Dict; hdutype=TableHDU, name=nothing, ver=nothing, header=nothing, units=nothing, varcols=nothing)\n\nCreate a new table extension and write data to it. If the FITS file is currently empty then a dummy primary array will be created before appending the table extension to it. data should be a dictionary with String keys (giving the column names) and Array values (giving data to write to each column). The following types are supported in binary tables: UInt8, Int8, UInt16, Int16, UInt32, Int32, Int64, Float32, Float64, Complex{Float32}, Complex{Float64}, String, Bool.\n\nOptional inputs:\n\nhdutype: Type of table extension to create. Can be either TableHDU (binary table) or ASCIITableHDU (ASCII table).\nname: Name of extension.\nver: Version of extension (Int).\nheader: FITSHeader instance to write to new extension.\nunits: Dictionary mapping column name to units (as a string).\nvarcols: An array giving the column names or column indicies to write as \"variable-length columns\".\n\nnote: Variable length columns\nVariable length columns allow a column\'s row entries to contain arrays of different lengths. They can potentially save diskspace when the rows of a column vary greatly in length, as the column data is all written to a contiguous heap area at the end of the table. Only column data of type Vector{String} or types such as Vector{Vector{UInt8}} can be written as variable length columns. In the second case, ensure that the column data type is a leaf type. That is, the type cannot be Vector{Vector{T}}, which would be an array of arrays having potentially non-uniform element types (which would not be writable as a FITS table column).\n\n\n\n\n\n"
},

{
    "location": "api.html#Base.write-Tuple{FITS,Array{String,1},Array{T,1} where T}",
    "page": "API Reference",
    "title": "Base.write",
    "category": "method",
    "text": "write(f::FITS, colnames, coldata; hdutype=TableHDU, name=nothing, ver=nothing, header=nothing, units=nothing, varcols=nothing)\n\nSame as write(f::FITS, data::Dict; ...) but providing column names and column data as a separate arrays. This is useful for specifying the order of the columns. Column names must be Vector{String} and column data must be a vector of arrays.\n\n\n\n\n\n"
},

{
    "location": "api.html#Base.read-Tuple{TableHDU,String}",
    "page": "API Reference",
    "title": "Base.read",
    "category": "method",
    "text": "read(hdu, colname; case_sensitive=true)\n\nRead a column as an array from the given table HDU.\n\nThe column name may contain wild card characters (*, ?, or #). The * wild card character matches any sequence of characters (including zero characters) and the ? character matches any single character. The # wildcard will match any consecutive string of decimal digits (0-9). The string must match a unique column.  The optional boolean keyword case_sensitive, true by default, specifies whether the column name is to be considered case sensitive.\n\n\n\n\n\n"
},

{
    "location": "api.html#Table-operations-1",
    "page": "API Reference",
    "title": "Table operations",
    "category": "section",
    "text": "FITSIO.colnames\nwrite(::FITS, ::Dict{String})\nwrite(::FITS, ::Vector{String}, ::Vector)\nread(::TableHDU, ::String)"
},

{
    "location": "api.html#FITSIO.libcfitsio_version",
    "page": "API Reference",
    "title": "FITSIO.libcfitsio_version",
    "category": "function",
    "text": "libcfitsio_version() -> VersionNumber\n\nReturn the version of the underlying CFITSIO library\n\nExample\n\njulia> FITSIO.libcfitsio_version()\nv\"3.37.0\"\n\n\n\n\n\n"
},

{
    "location": "api.html#Miscellaneous-1",
    "page": "API Reference",
    "title": "Miscellaneous",
    "category": "section",
    "text": "FITSIO.libcfitsio_version"
},

{
    "location": "libcfitsio.html#",
    "page": "Libcfitsio Submodule",
    "title": "Libcfitsio Submodule",
    "category": "page",
    "text": ""
},

{
    "location": "libcfitsio.html#Libcfitsio-submodule-1",
    "page": "Libcfitsio Submodule",
    "title": "Libcfitsio submodule",
    "category": "section",
    "text": "CurrentModule = FITSIO.LibcfitsioThe Libcfitsio submodule provides an interface familiar to users of the CFITSIO C library. It can be used withusing FITSIO.LibcfitsioThe functions exported by this module operate on FITSFile objects, which is a thin wrapper around a pointer to a CFITSIO fitsfile.  For the most part, the functions are thin wrappers around the CFITSIO routines of the same names. Typically, they:Convert from Julia types to C types as necessary.\nCheck the returned status value and raise an appropriate exception if non-zero.warning: Warning\nNote that these functions do not check if the file is still open before trying to access it. A segmentation fault can result from trying to operate on a closed file. (The main FITSIO interface always checks if the file is open before any operation.)"
},

{
    "location": "libcfitsio.html#FITSIO.Libcfitsio.fits_create_file",
    "page": "Libcfitsio Submodule",
    "title": "FITSIO.Libcfitsio.fits_create_file",
    "category": "function",
    "text": "fits_create_file(filename::AbstractString)\n\nCreate and open a new empty output FITSFile.\n\n\n\n\n\n"
},

{
    "location": "libcfitsio.html#FITSIO.Libcfitsio.fits_clobber_file",
    "page": "Libcfitsio Submodule",
    "title": "FITSIO.Libcfitsio.fits_clobber_file",
    "category": "function",
    "text": "fits_clobber_file(filename::AbstractString)\n\nLike fits_create_file, but overwrites filename if it exists.\n\n\n\n\n\n"
},

{
    "location": "libcfitsio.html#FITSIO.Libcfitsio.fits_open_file",
    "page": "Libcfitsio Submodule",
    "title": "FITSIO.Libcfitsio.fits_open_file",
    "category": "function",
    "text": "fits_open_file(filename::String)\n\nOpen an existing data file.\n\n\n\n\n\n"
},

{
    "location": "libcfitsio.html#FITSIO.Libcfitsio.fits_open_table",
    "page": "Libcfitsio Submodule",
    "title": "FITSIO.Libcfitsio.fits_open_table",
    "category": "function",
    "text": "fits_open_table(filename::String)\n\nOpen an existing data file (like fits_open_file) and move to the first HDU containing either an ASCII or a binary table.\n\n\n\n\n\n"
},

{
    "location": "libcfitsio.html#FITSIO.Libcfitsio.fits_open_image",
    "page": "Libcfitsio Submodule",
    "title": "FITSIO.Libcfitsio.fits_open_image",
    "category": "function",
    "text": "fits_open_image(filename::String)\n\nOpen an existing data file (like fits_open_file) and move to the first HDU containing an image.\n\n\n\n\n\n"
},

{
    "location": "libcfitsio.html#FITSIO.Libcfitsio.fits_open_data",
    "page": "Libcfitsio Submodule",
    "title": "FITSIO.Libcfitsio.fits_open_data",
    "category": "function",
    "text": "fits_open_data(filename::String)\n\nOpen an existing data file (like fits_open_file) and move to the first HDU containing either an image or a table.\n\n\n\n\n\n"
},

{
    "location": "libcfitsio.html#FITSIO.Libcfitsio.fits_close_file",
    "page": "Libcfitsio Submodule",
    "title": "FITSIO.Libcfitsio.fits_close_file",
    "category": "function",
    "text": "fits_close_file(f::FITSFile)\n\nClose a previously opened FITS file.\n\n\n\n\n\n"
},

{
    "location": "libcfitsio.html#FITSIO.Libcfitsio.fits_delete_file",
    "page": "Libcfitsio Submodule",
    "title": "FITSIO.Libcfitsio.fits_delete_file",
    "category": "function",
    "text": "fits_delete_file(f::FITSFile)\n\nClose an opened FITS file (like fits_close_file) and removes it from the disk.\n\n\n\n\n\n"
},

{
    "location": "libcfitsio.html#FITSIO.Libcfitsio.fits_file_name",
    "page": "Libcfitsio Submodule",
    "title": "FITSIO.Libcfitsio.fits_file_name",
    "category": "function",
    "text": "fits_file_name(f::FITSFile)\n\nReturn the name of the file associated with object f.\n\n\n\n\n\n"
},

{
    "location": "libcfitsio.html#File-access-1",
    "page": "Libcfitsio Submodule",
    "title": "File access",
    "category": "section",
    "text": "fits_create_file\nfits_clobber_file\nfits_open_file\nfits_open_table\nfits_open_image\nfits_open_data\nfits_close_file\nfits_delete_file\nfits_file_name"
},

{
    "location": "libcfitsio.html#FITSIO.Libcfitsio.fits_get_num_hdus",
    "page": "Libcfitsio Submodule",
    "title": "FITSIO.Libcfitsio.fits_get_num_hdus",
    "category": "function",
    "text": "fits_get_num_hdus(f::FITSFile)\n\nReturn the number of HDUs in the file.\n\n\n\n\n\n"
},

{
    "location": "libcfitsio.html#FITSIO.Libcfitsio.fits_movabs_hdu",
    "page": "Libcfitsio Submodule",
    "title": "FITSIO.Libcfitsio.fits_movabs_hdu",
    "category": "function",
    "text": "fits_movabs_hdu(f::FITSFile, hduNum::Integer)\n\nChange the current HDU to the value specified by hduNum, and return a symbol describing the type of the HDU.\n\nPossible symbols are: image_hdu, ascii_table, or binary_table. The value of hduNum must range between 1 and the value returned by fits_get_num_hdus.\n\n\n\n\n\n"
},

{
    "location": "libcfitsio.html#FITSIO.Libcfitsio.fits_movrel_hdu",
    "page": "Libcfitsio Submodule",
    "title": "FITSIO.Libcfitsio.fits_movrel_hdu",
    "category": "function",
    "text": "fits_movrel_hdu(f::FITSFile, hduNum::Integer)\n\nChange the current HDU by moving forward or backward by hduNum HDUs (positive means forward), and return the same as fits_movabs_hdu.\n\n\n\n\n\n"
},

{
    "location": "libcfitsio.html#FITSIO.Libcfitsio.fits_movnam_hdu",
    "page": "Libcfitsio Submodule",
    "title": "FITSIO.Libcfitsio.fits_movnam_hdu",
    "category": "function",
    "text": "fits_movnam_hdu(f::FITSFile, extname::String, extver::Integer=0,\n                hdu_type_int::Integer=-1)\n\nChange the current HDU by moving to the (first) HDU which has the specified extension type and EXTNAME and EXTVER keyword values (or HDUNAME and HDUVER keywords).\n\nIf extver is 0 (the default) then the EXTVER keyword is ignored and the first HDU with a matching EXTNAME (or HDUNAME) keyword will be found. If hdu_type_int is -1 (the default) only the extname and extver values will be used to locate the correct extension. If no matching HDU is found in the file, the current HDU will remain unchanged.\n\n\n\n\n\n"
},

{
    "location": "libcfitsio.html#HDU-Routines-1",
    "page": "Libcfitsio Submodule",
    "title": "HDU Routines",
    "category": "section",
    "text": "The functions described in this section change the current HDU and to find their number and type. The following is a short example which shows how to use them:num = fits_get_num_hdus(f)\nprintln(\"Number of HDUs in the file: \", num)\n\nfor i = 1:num\n    hdu_type = fits_movabs_hdu(f, i)\n    println(i, \") hdu_type = \", hdu_type)\nendfits_get_num_hdus\nfits_movabs_hdu\nfits_movrel_hdu\nfits_movnam_hdu"
},

{
    "location": "libcfitsio.html#FITSIO.Libcfitsio.fits_get_hdrspace",
    "page": "Libcfitsio Submodule",
    "title": "FITSIO.Libcfitsio.fits_get_hdrspace",
    "category": "function",
    "text": "fits_get_hdrspace(f::FITSFile) -> (keysexist, morekeys)\n\nReturn the number of existing keywords (not counting the END keyword) and the amount of space currently available for more keywords.\n\n\n\n\n\n"
},

{
    "location": "libcfitsio.html#FITSIO.Libcfitsio.fits_read_keyword",
    "page": "Libcfitsio Submodule",
    "title": "FITSIO.Libcfitsio.fits_read_keyword",
    "category": "function",
    "text": "fits_read_keyword(f::FITSFile, keyname::String) -> (value, comment)\n\nyields the specified keyword value and commend (as a tuple of strings), throws and error if the keyword is not found.\n\n\n\n\n\n"
},

{
    "location": "libcfitsio.html#FITSIO.Libcfitsio.fits_read_record",
    "page": "Libcfitsio Submodule",
    "title": "FITSIO.Libcfitsio.fits_read_record",
    "category": "function",
    "text": "fits_read_record(f::FITSFile, keynum::Int) -> String\n\nReturn the nth header record in the CHU. The first keyword in the header is at keynum = 1.\n\n\n\n\n\n"
},

{
    "location": "libcfitsio.html#FITSIO.Libcfitsio.fits_read_keyn",
    "page": "Libcfitsio Submodule",
    "title": "FITSIO.Libcfitsio.fits_read_keyn",
    "category": "function",
    "text": "fits_read_keyn(f::FITSFile, keynum::Int) -> (name, value, comment)\n\nReturn the nth header record in the CHU. The first keyword in the header is at keynum = 1.\n\n\n\n\n\n"
},

{
    "location": "libcfitsio.html#FITSIO.Libcfitsio.fits_write_key",
    "page": "Libcfitsio Submodule",
    "title": "FITSIO.Libcfitsio.fits_write_key",
    "category": "function",
    "text": "fits_write_key(f::FITSFile, keyname::String, value, comment::String)\n\nWrite a keyword of the appropriate data type into the CHU.\n\n\n\n\n\n"
},

{
    "location": "libcfitsio.html#FITSIO.Libcfitsio.fits_write_record",
    "page": "Libcfitsio Submodule",
    "title": "FITSIO.Libcfitsio.fits_write_record",
    "category": "function",
    "text": "fits_write_record(f::FITSFile, card::String)\n\nWrite a user specified keyword record into the CHU.\n\n\n\n\n\n"
},

{
    "location": "libcfitsio.html#FITSIO.Libcfitsio.fits_delete_record",
    "page": "Libcfitsio Submodule",
    "title": "FITSIO.Libcfitsio.fits_delete_record",
    "category": "function",
    "text": "fits_delete_record(f::FITSFile, keynum::Int)\n\nDelete the keyword record at the specified index.\n\n\n\n\n\n"
},

{
    "location": "libcfitsio.html#FITSIO.Libcfitsio.fits_delete_key",
    "page": "Libcfitsio Submodule",
    "title": "FITSIO.Libcfitsio.fits_delete_key",
    "category": "function",
    "text": "fits_delete_key(f::FITSFile, keyname::String)\n\nDelete the keyword named keyname.\n\n\n\n\n\n"
},

{
    "location": "libcfitsio.html#FITSIO.Libcfitsio.fits_hdr2str",
    "page": "Libcfitsio Submodule",
    "title": "FITSIO.Libcfitsio.fits_hdr2str",
    "category": "function",
    "text": "fits_hdr2str(f::FITSFile, nocomments::Bool=false)\n\nReturn the header of the CHDU as a string. If nocomments is true, comment cards are stripped from the output.\n\n\n\n\n\n"
},

{
    "location": "libcfitsio.html#Header-Keyword-Routines-1",
    "page": "Libcfitsio Submodule",
    "title": "Header Keyword Routines",
    "category": "section",
    "text": "fits_get_hdrspace\nfits_read_keyword\nfits_read_record\nfits_read_keyn\nfits_write_key\nfits_write_record\nfits_delete_record\nfits_delete_key\nfits_hdr2str"
},

{
    "location": "libcfitsio.html#FITSIO.Libcfitsio.fits_get_img_size",
    "page": "Libcfitsio Submodule",
    "title": "FITSIO.Libcfitsio.fits_get_img_size",
    "category": "function",
    "text": "fits_get_img_size(f::FITSFile)\n\nGet the dimensions of the image.\n\n\n\n\n\n"
},

{
    "location": "libcfitsio.html#FITSIO.Libcfitsio.fits_create_img",
    "page": "Libcfitsio Submodule",
    "title": "FITSIO.Libcfitsio.fits_create_img",
    "category": "function",
    "text": "fits_create_img(f::FITSFile, t::Type, naxes::Vector{Int})\n\nCreate a new primary array or IMAGE extension with a specified data type and size.\n\n\n\n\n\n"
},

{
    "location": "libcfitsio.html#FITSIO.Libcfitsio.fits_write_pix",
    "page": "Libcfitsio Submodule",
    "title": "FITSIO.Libcfitsio.fits_write_pix",
    "category": "function",
    "text": "fits_write_pix(f::FITSFile, fpixel::Vector{Int}, nelements::Int, data::Array)\n\nWrite pixels from data into the FITS file.\n\n\n\n\n\n"
},

{
    "location": "libcfitsio.html#FITSIO.Libcfitsio.fits_read_pix",
    "page": "Libcfitsio Submodule",
    "title": "FITSIO.Libcfitsio.fits_read_pix",
    "category": "function",
    "text": "fits_read_pix(f::FITSFile, fpixel::Vector{Int}, nelements::Int, data::Array)\n\nRead pixels from the FITS file into data.\n\n\n\n\n\n"
},

{
    "location": "libcfitsio.html#Image-HDU-Routines-1",
    "page": "Libcfitsio Submodule",
    "title": "Image HDU Routines",
    "category": "section",
    "text": "fits_get_img_size\nfits_create_img\nfits_write_pix\nfits_read_pix"
},

{
    "location": "libcfitsio.html#FITSIO.Libcfitsio.fits_create_ascii_tbl",
    "page": "Libcfitsio Submodule",
    "title": "FITSIO.Libcfitsio.fits_create_ascii_tbl",
    "category": "function",
    "text": "fits_create_ascii_tbl(f::FITSFile, numrows::Integer, coldefs::Array{ColumnDef},\n                      extname::String)\n\nAppend a new HDU containing an ASCII table.\n\nThe table will have numrows rows (this parameter can be set to zero), each initialized with the default value. In order to create a table, the programmer must specify the characteristics of each column. The columns are specified by the coldefs variable, which is an array of tuples. Each tuple must have three string fields:\n\nThe name of the column.\nThe data type and the repetition count. It must be a string made by a number (the repetition count) followed by a letter specifying the type (in the example above, D stands for Float64, E stands for Float32, A stands for Char). Refer to the CFITSIO documentation for more information about the syntax of this parameter.\nThe measure unit of this field. This is used only as a comment.\n\nThe value of extname sets the \"extended name\" of the table, i.e., a string that in some situations can be used to refer to the HDU itself.\n\nNote that, unlike for binary tables, CFITSIO puts some limitations to the types that can be used in an ASCII table column. Refer to the CFITSIO manual for further information.\n\nSee also fits_create_binary_tbl for a similar function which creates binary tables.\n\n\n\n\n\n"
},

{
    "location": "libcfitsio.html#FITSIO.Libcfitsio.fits_create_binary_tbl",
    "page": "Libcfitsio Submodule",
    "title": "FITSIO.Libcfitsio.fits_create_binary_tbl",
    "category": "function",
    "text": "fits_create_binary_tbl(f::FITSFile, numrows::Integer, coldefs::Array{ColumnDef},\n                       extname::String)\n\nAppend a new HDU containing a binary table. The meaning of the parameters is the same as in a call to fits_create_ascii_tbl.\n\nIn general, one should pick this function for creating tables in a new HDU, as binary tables require less space on the disk and are more efficient to read and write. (Moreover, a few datatypes are not supported in ASCII tables).\n\n\n\n\n\n"
},

{
    "location": "libcfitsio.html#FITSIO.Libcfitsio.fits_get_coltype",
    "page": "Libcfitsio Submodule",
    "title": "FITSIO.Libcfitsio.fits_get_coltype",
    "category": "function",
    "text": "fits_get_coltype(f::FITSFile, colnum::Integer)\n\nProvided that the current HDU contains either an ASCII or binary table, return information about the column at position colnum (counting from 1).\n\nReturn is a tuple containing\n\ntypecode: CFITSIO integer type code of the column.\nrepcount: Repetition count for the column.\nwidth: Width of an individual element.\n\n\n\n\n\n"
},

{
    "location": "libcfitsio.html#FITSIO.Libcfitsio.fits_insert_rows",
    "page": "Libcfitsio Submodule",
    "title": "FITSIO.Libcfitsio.fits_insert_rows",
    "category": "function",
    "text": "fits_insert_rows(f::FITSFile, firstrow::Integer, nrows::Integer)\n\nInsert a number of rows equal to nrows after the row number firstrow.\n\nThe elements in each row are initialized to their default value: you can modify them later using fits_write_col.\n\nSince the first row is at position 1, in order to insert rows before the first one firstrow must be equal to zero.\n\n\n\n\n\n"
},

{
    "location": "libcfitsio.html#FITSIO.Libcfitsio.fits_delete_rows",
    "page": "Libcfitsio Submodule",
    "title": "FITSIO.Libcfitsio.fits_delete_rows",
    "category": "function",
    "text": "fits_delete_rows(f::FITSFile, firstrow::integer, nrows::Integer)\n\nDelete nrows rows, starting from the one at position firstrow. The index of the first row is 1.\n\n\n\n\n\n"
},

{
    "location": "libcfitsio.html#FITSIO.Libcfitsio.fits_read_col",
    "page": "Libcfitsio Submodule",
    "title": "FITSIO.Libcfitsio.fits_read_col",
    "category": "function",
    "text": "fits_read_col(f, colnum, firstrow, firstelem, data)\n\nRead data from one column of an ASCII/binary table and convert the data into the specified type T.\n\nArguments\n\nf::FITSFile: the file to be read.\ncolnum::Integer: the column number, where the value of the first column is 1.\nfirstrow::Integer: the elements to be read start from this row.\nfirstelem::Integer: specifies which is the first element to be read, when each cell contains more than one element (i.e., the \"repetition count\" of the field is greater than one).\ndata::Array: at the end of the call, this will be filled with the elements read from the column. The length of the array gives the overall number of elements.\n\n\n\n\n\n"
},

{
    "location": "libcfitsio.html#FITSIO.Libcfitsio.fits_write_col",
    "page": "Libcfitsio Submodule",
    "title": "FITSIO.Libcfitsio.fits_write_col",
    "category": "function",
    "text": "fits_write_col(f, colnum, firstrow, firstelem, data)\n\nWrite some data in one column of a ASCII/binary table.\n\nIf there is no room for the elements, new rows will be created. (It is therefore useless to call fits_insert_rows if you only need to append elements to the end of a table.)\n\nf::FITSFile: the file in which data will be written.\ncolnum::Integer: the column number, where the value of the first column is 1.\nfirstrow::Integer: the data wil be written from this row onwards.\nfirstelem::Integer: specifies the position in the row where the first element will be written.\ndata::Array: contains the elements that are to be written to the column of the table.\n\n\n\n\n\n"
},

{
    "location": "libcfitsio.html#Table-Routines-1",
    "page": "Libcfitsio Submodule",
    "title": "Table Routines",
    "category": "section",
    "text": "There are two functions to create a new HDU table extension: fits_create_ascii_table and fits_create_binary_table. In general, one should pick the second as binary tables require less space on the disk and are more efficient to read and write. (Moreover, a few datatypes are not supported in ASCII tables). In order to create a table, the programmer must specify the characteristics of each column by passing an array of tuples. Here is an example:f = fits_create_file(\"!new.fits\")\ncoldefs = [(\"SPEED\", \"1D\", \"m/s\"),\n           (\"MASS\", \"1E\", \"kg\"),\n           (\"PARTICLE\", \"20A\", \"Name\")]\nfits_create_binary_tbl(f, 10, coldefs, \"PARTICLE\")This example creates a table with room for 10 entries, each of them describing the characteristics of a particle: its speed, its mass, and its name (codified as a 20-character string). See the documentation of fits_create_ascii_tbl for more details.fits_create_ascii_tbl\nfits_create_binary_tbl\nfits_get_coltype\nfits_insert_rows\nfits_delete_rows\nfits_read_col\nfits_write_col"
},

]}
