module FITSIO

export FITSFile,
       fits_clobber_file,
       fits_close_file,
       fits_copy_image_section,
       fits_create_ascii_tbl,
       fits_create_binary_tbl,
       fits_create_file,
       fits_create_img,
       fits_delete_file,
       fits_delete_key,
       fits_delete_record,
       fits_delete_rows,
       fits_file_mode,
       fits_file_name,
       fits_get_col_repeat,
       fits_get_hdrspace,
       fits_get_hdu_num,
       fits_get_hdu_type,
       fits_get_img_dim,
       fits_get_img_equivtype,
       fits_get_img_size,
       fits_get_img_type,
       fits_get_num_cols,
       fits_get_num_hdus,
       fits_get_num_rows,
       fits_get_num_rowsll,
       fits_get_rowsize,
       fits_hdr2str,
       fits_insert_rows,
       fits_movabs_hdu,
       fits_movrel_hdu,
       fits_movnam_hdu,
       fits_open_data,
       fits_open_file,
       fits_open_image,
       fits_open_table,
       fits_read_col,
       fits_read_keyn,
       fits_read_keyword,
       fits_read_pix,
       fits_read_record,
       fits_read_subset,
       fits_write_col,
       fits_write_key,
       fits_write_pix,
       fits_write_record

# HDU type interface
export FITS,
       HDU,
       ImageHDU,
       FITSHeader,
       readkey,
       readheader,
       getcomment,
       setcomment!,
       copy_section

import Base: getindex, setindex!, length, show, read, write, close, ndims,
             size, endof, haskey, keys

using BinDeps
@BinDeps.load_dependencies

include("cfitsio.jl")  # Low-level cfitsio functions
include("hdutypes.jl")  # HDU type interface

end # module
