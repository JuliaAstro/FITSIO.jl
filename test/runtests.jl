using Aqua
using CFITSIO
using FITSIO
using OrderedCollections
import Tables

# Deal with compatibility issues.
using Test
using Random # for `randstring`

@testset "project quality" begin
    Aqua.test_all(FITSIO)
end

randomfilename(dir) = joinpath(dir, randstring(12)*".fits")
function tempnamefits(f)
    mktempdir() do dir
        f(randomfilename(dir))
    end
end

@testset "File IO" begin
    tempnamefits() do fname
        # create the file first
        FITS(fname, "w") do f
            write(f, ones(2))
        end
        # open as read only
        f = FITS(fname, "r")
        # writing should throw an error
        @test_throws Exception write(f, ones(2))
        d = f[1]
        @test_throws Exception write(d, ones(2))

        fname2 = fname*"[12].fits"
        FITS(fname2, "w", extendedparser = false) do f
            write(f, [1,2])
            @test read(f[1]) == [1,2]
        end
        FITSIO.fitswrite(fname2, [1:4;], extendedparser = false)
        @test FITSIO.fitsread(fname2, extendedparser = false) == [1:4;]
    end
end

@testset "Images" begin
    # Create a FITS instance and loop over supported types.
    tempnamefits() do fname
        FITS(fname, "w") do f
            for T in (UInt8, Int8, UInt16, Int16, UInt32, Int32, Int64,
                      Float32, Float64)
                indata = reshape(T[1:100;], 5, 20)

                # Test writing the data to a new extension
                write(f, indata)

                # test reading the full array
                outdata = @inferred read(f[end])
                @test indata == outdata
                @test eltype(indata) == eltype(outdata) == eltype(f[end])

                # test reading subsets of the array
                @inferred read(f[end], :, :) == indata
                @inferred read(f[end], 4, 1:10) == indata[4, 1:10]  # 2-d array
                @inferred read(f[end], :, 4) == indata[:, 4]  # 1-d array
                @inferred read(f[end], 2, 3) == indata[2, 3]  # scalar
                @inferred read(f[end], :, 1:2:10) == indata[:, 1:2:10]
                @inferred read(f[end], 1:3, :) == indata[1:3, :]

                # test expected errors
                @test_throws DimensionMismatch read(f[end], :)
                @test_throws DimensionMismatch read(f[end], :, :, 1)
                @test_throws BoundsError read(f[end], 1:6, :)
                @test_throws BoundsError read(f[end], 1, 0)

            end

            @test_throws ErrorException f[100]

            # Test representation
            @test repr(f)[end-17:end] == "9          Image  "
            @test repr(f[1])[1:6] == "File: "

            # test iteration
            for hdu in f
                @inferred size(hdu) == (5, 20)
            end
        end
    end

    @testset "copy_section" begin
        mktempdir() do dir
            fname1 = randomfilename(dir)
            FITS(fname1, "w") do f1
                indata = reshape(Float32[1:400;], 20, 20)
                write(f1, indata)

                fname2 = randomfilename(dir)
                FITS(fname2, "w") do f2
                    copy_section(f1[1], f2, 1:10, 1:10)
                    copy_section(f1[1], f2, 1:10, 1:2:20)
                    outdata = read(f2[1])
                    @test outdata == indata[1:10, 1:10]
                    outdata = read(f2[2])
                    @test outdata == indata[1:10, 1:2:20]
                end
            end
        end
    end

    @testset "non-allocating read" begin
        tempnamefits() do fname
            FITS(fname, "r+") do f

                indata = reshape(Float32[1:400;], 20, 20)
                write(f, indata)

                indata3D = reshape([1:30;], 3, 2, 5)
                write(f, indata3D)

                # Read the entire array using read to compare with read!
                a = read(f[1])
                a3D = read(f[2])

                # Read the entire array
                b = zeros(eltype(indata),size(indata))
                read!(f[1],b)
                @test a == b
                read!(f[1],b,:,:)
                @test a == b

                # Read a 2D slice
                b = zeros(eltype(indata),2,2)
                read!(f[1],b,1:2,1:2)
                @test a[1:2,1:2] == b

                # Read an entire 1D slice
                b = zeros(eltype(indata),size(indata,1))
                read!(f[1],b,:,1)
                @test a[:,1] == b

                b = zeros(eltype(indata),size(indata,2))
                read!(f[1],b,1,:)
                @test a[1,:] == b

                # Read a part of a 1D slice
                b = zeros(eltype(indata),1)
                read!(f[1],b,1:1,1)
                @test a[1,1] == b[1]
                read!(f[1],b,1,1:1)
                @test a[1,1] == b[1]

                # Read a single element into a 0-dim array
                b = zeros(eltype(indata))
                read!(f[1],b,1,1)
                @test a[1,1] == first(b)
                read!(f[2], b, 1, 1, 1)
                @test a3D[1,1,1] == first(b)
                # read the entire image into an array with a different eltype
                b3D = zeros(size(a3D))
                read!(f[2], b3D)
                @test a3D == b3D

                b = zeros(eltype(indata))
                @test_throws DimensionMismatch read!(f[1],b,1:10,1)
                @test_throws DimensionMismatch read!(f[1],b,1:10,1:10)
                @test_throws DimensionMismatch read!(f[1],b)
                @test_throws DimensionMismatch read!(f[1],b,:,1)
                @test_throws DimensionMismatch read!(f[1],b,1,:)
                @test_throws DimensionMismatch read!(f[1],b,:,:)
                @test_throws DimensionMismatch read!(f[1], zeros(1,1))
                @test_throws DimensionMismatch read!(f[1], zeros(1,1), :)

                b3D = zero(indata3D)
                read!(f[2],b3D)
                @test a3D == b3D
                read!(f[2],b3D,:,:,:)
                @test a3D == b3D

                b0D = zeros(eltype(indata))
                read!(f[1],b0D,1,1)
                @test a[1,1] == b0D[]

                @testset "different ndims" begin
                    local b = similar(indata, length(indata))
                    read!(f[1], b)
                    @test vec(b) == vec(a)
                end

                @testset "read into views" begin

                    local b = zero(indata)

                    # Entire array
                    b_view = @view b[:,:];
                    read!(f[1], b_view);
                    @test a == b

                    b_view = @view b3D[:,:,:];
                    read!(f[2], b_view);
                    @test a3D == b3D;

                    @testset "1D slices of a 2D array" begin
                        for ax2 in axes(b,2)
                            b .= zero(eltype(b))
                            b_view = @view b[:,ax2]
                            read!(f[1],b_view,:,ax2)
                            @test a[:,ax2] == b[:,ax2]
                        end

                        # Non-contiguous views can not be read into
                        b .= zero(eltype(b))
                        b_view = @view b[1,:]
                        @test_throws ArgumentError read!(f[1], b_view, 1, :)
                    end

                    @testset "1D slices of a 3D array" begin
                        for ax2 in axes(b3D,2), ax3 in axes(b3D,3)
                            b3D .= zero(eltype(b3D))
                            b_view = @view b3D[:,ax2,ax3]
                            read!(f[2],b_view,:,ax2,ax3)
                            @test a3D[:,ax2,ax3] == b3D[:,ax2,ax3]
                        end

                        # Non-contiguous views can not be read into
                        b3D .= zero(eltype(b3D))
                        b_view = @view b3D[1,1,:]
                        @test_throws ArgumentError read!(f[2],b_view,1,1,:)

                        b_view = @view b3D[1,:,1]
                        @test_throws ArgumentError read!(f[2],b_view,1,:,1)
                    end

                    @testset "2D slices of a 3D array" begin
                        for ax3 in axes(b3D,3)
                            b3D .= zero(eltype(b3D))
                            b_view = @view b3D[:,:,ax3]
                            read!(f[2],b_view,:,:,ax3)
                            @test a3D[:,:,ax3] == b3D[:,:,ax3]
                        end

                        # Non-contiguous views can not be read into
                        b3D .= zero(eltype(b3D))
                        b_view = @view b3D[:,1,:]
                        @test_throws ArgumentError read!(f[2],b_view,:,1,:)

                        b_view = @view b3D[1,:,:]
                        @test_throws ArgumentError read!(f[2],b_view,1,:,:)
                    end

                    @testset "0D view" begin
                        b .= zero(eltype(b))
                        b_view = @view b[1,1]
                        read!(f[1],b_view,1,1)
                        @test a[1,1] == b[1,1]
                    end
                end
            end
        end
    end

    @testset "reinterpreted complex" begin
        @testset "read" begin
            function _testreadcomplex(fname, arr::AbstractArray{Complex{T}}) where T
                FITS(fname,"w") do f
                    write(f, reinterpret(T,arr))
                end

                b = similar(arr)
                FITS(fname,"r") do f
                    read!(f[1], reinterpret(T,b))
                end
                @test b == arr

                c = FITS(fname,"r") do f
                    read(f[1])
                end
                @test c == reinterpret(T,arr)
            end

            function testreadcomplex(fname, arr)
                # Given a complex array, reinterpret it as float and write to file
                # Read it back and compare
                _testreadcomplex(fname, arr)

                # Write a contiguous subsection instead of the entire array
                region = 1:size(arr,1)
                arrv = @view arr[region]
                _testreadcomplex(fname, arrv)
            end

            tempnamefits() do fname
                testreadcomplex(fname, ones(ComplexF64,3))
                testreadcomplex(fname, ones(ComplexF64,3,4))
                testreadcomplex(fname, ones(ComplexF64,3,4,5))
            end
        end

        @testset "write" begin
            function _testwritecomplex(fname, arr::AbstractArray{Complex{T}}) where T
                FITS(fname,"w") do f
                    write(f, reinterpret(T,arr) )
                end

                a = FITS(fname,"r") do f
                    read(f[1])
                end

                b = reinterpret(T,arr)
                @test b == a

                c = FITS(fname,"r") do f
                    reinterpret(Complex{T}, read(f[1]))
                end
                @test c == arr
            end

            function testwritecomplex(fname, arr)
                # Given a complex array, reinterpret it as float and write to file
                # Read it back and compare
                _testwritecomplex(fname, arr)

                # Write a subsection instead of the entire array
                region = 1:size(arr,1)
                arrv = @view arr[region]
                _testwritecomplex(fname, arrv)
            end

            tempnamefits() do fname
                testwritecomplex(fname, ones(ComplexF64,3))
                testwritecomplex(fname, ones(ComplexF64,3,4))
                testwritecomplex(fname, ones(ComplexF64,3,4,5))
            end
        end
    end

    @testset "write views to file" begin
        b = reshape([1:8;], 2, 2, 2)

        @testset "2D slice of a 3D array" begin
            tempnamefits() do fname
                FITS(fname, "r+") do f
                    for (ind,ax3) in enumerate(axes(b,3))
                        b_view = @view b[:,:,ax3]
                        write(f,b_view)
                        @test read(f[ind]) == b_view
                    end
                end
            end
        end

        @testset "1D slice of a 3D array" begin
            tempnamefits() do fname
                FITS(fname, "r+") do f
                    ax23iter = Iterators.product(axes(b)[2:3]...)
                    for (ind,(ax2,ax3)) in enumerate(ax23iter)
                        b_view = @view b[:,ax2,ax3]
                        write(f,b_view)
                        @test read(f[ind]) == b_view
                    end
                end
            end
        end

        @testset "Non-contiguous" begin
            tempnamefits() do fname
                FITS(fname, "r+") do f
                    b_view = @view b[:,1,:]
                    @test_throws ArgumentError write(f,b_view)

                    b_view = @view b[1,:,:]
                    @test_throws ArgumentError write(f,b_view)

                    b_view = @view b[1,1,:]
                    @test_throws ArgumentError write(f,b_view)

                    b_view = @view b[1,:,1]
                    @test_throws ArgumentError write(f,b_view)

                    # write something random to avoid the error on closing
                    write(f,zeros(1))
                end
            end
        end
    end

    @testset "fitsread" begin
        a = ones(3,3)

        tempnamefits() do fname
            FITS(fname, "w") do f
                write(f, a)
            end
            @test FITSIO.fitsread(fname) == a
            @test FITSIO.fitsread(fname, 1) == a
        end
    end

    # Ref: https://github.com/JuliaAstro/FITSIO.jl/issues/163
    @testset "String key" begin
        tempnamefits() do fname
            # Create test data
            a = ones(3,3)
            
            FITS(fname, "w") do f
                GC.@preserve a begin
                    # Write data with name
                    write(f, a, name="a")
                    
                    # For Julia < 1.6, handle strings explicitly
                    @static if VERSION < v"1.6"
                        key_a = String("a")
                        GC.@preserve key_a begin
                            read_data = read(f[key_a])
                            @test read_data == a
                        end
                        
                        # Read by index
                        @test read(f[1]) == a
                        
                        # Test haskey with strings
                        key_a = String("a")
                        key_b = String("b")
                        GC.@preserve key_a key_b begin
                            @test haskey(f, key_a)
                            @test !haskey(f, key_b)
                        end
                    else
                        @test read(f["a"]) == a
                        @test read(f[1]) == a
                        @test haskey(f, "a")
                        @test !haskey(f, "b")
                    end
                    
                    # Test numeric keys (same for all versions)
                    @test haskey(f, 1)
                    @test !haskey(f, 2)
                end
            end
        end
    end

    @testset "non-Int integer indices" begin
        tempnamefits() do filename
            FITS(filename, "w") do f
                write(f, ones(2,2))
                @test read(f[1], big(1), Int8(1)) == read(f[1], 1, 1)
                @test read(f[1], big(1), Int8(1):Int8(2)) == read(f[1], 1, 1:2)
                @test read(f[1], :, Int8(1)) == read(f[1], :, 1)
            end
        end
    end

    @testset "delete" begin
        tempnamefits() do fname
            FITS(fname, "w") do f
                write(f, ones(2,2))
                write(f, ones(3,3))
                @test length(f) == 2
                hdu = f[2]
                deleteat!(f, 2)
                @test length(f) == 1
                @test repr(hdu) == "Deleted HDU"

                # if the array is read in before deletion,
                # test that the hdu is removed from the cache
                write(f, ones(3,3))
                read(f[2])
                deleteat!(f, 2)
                @test length(f) == 1

                write(f, ones(3,3))
                deleteat!(f,1)
                @test length(f) == 2
                @test ndims(f[1]) == 0
                write(f, ones(4,4))
                deleteat!(f, 2)
                @test size(f[2]) == (4,4)
            end
        end
    end
end

@testset "Write data to an existing image HDU" begin
    @testset "Overwrite an image" begin
        # Create a file with two images
        tempnamefits() do fname

            FITS(fname, "w") do f
                write(f, [[1 2 3]; [4 5 6]])
                write(f, [[7 8 9]; [17 52 10]])
            end

            # Open the file in read-write mode
            FITS(fname, "r+") do f
                @test read(f[1]) == [[1 2 3]; [4 5 6]]
                @test read(f[2]) == [[7 8 9]; [17 52 10]]

                # Write data into the first image HDU
                write(f[1], [[11 12 13]; [14 15 16]])

                @test read(f[1]) == [[11 12 13]; [14 15 16]]  # First HDU is updated
                @test read(f[2]) == [[7 8 9]; [17 52 10]]  # Second HDU is untouched
            end
        end
    end

    @testset "Show size mismatch error" begin
        tempnamefits() do fname
            FITS(fname, "w") do f
                write(f, [[1 2 3 4]; [5 6 7 8]])
                write(f, [[1 2 3]; [4 5 6]])
            end

            FITS(fname, "r+") do f
                second_hdu = f[2]  # Get the second HDU

                # Make the first HDU to be the current HDU
                # Needed to check if `write` method compares data size with the given HDU
                # and not with the current HDU
                fits_movabs_hdu(f.fitsfile, 1)

                # Attempt to write data array of incorrect size to the second HDU
                exception = ErrorException("size of HDU (2, 3) not equal to size of data (2, 4).")
                @test_throws exception write(second_hdu, [[11 12 13 14]; [15 16 17 18]])
            end
        end
    end

    @testset "write with different eltype" begin
        tempnamefits() do fname
            FITS(fname, "w") do f
                write(f, ones(2,2)) # create an image with eltype Float64
                write(f[1], ones(Int, 2, 2) .* 2) # write Int data to the same HDU
                @test read(f[1]) == ones(2,2) .* 2 # check that data is written correctly
            end
        end
    end

    @testset "fitswrite" begin
        tempnamefits() do fname
            a = ones(3,3)
            FITSIO.fitswrite(fname, a)
            FITS(fname, "r") do f
                @test read(f[1]) == a
            end
        end
    end
end

@testset "Tables" begin
    tempnamefits() do fname
        FITS(fname, "w") do f
            ## Binary table
            indata = Dict{String, Array}()
            for (i, T) in enumerate([UInt8, Int8, UInt16, Int16, UInt32, Int32, UInt64, Int64,
                                     Float32, Float64, ComplexF32, ComplexF64])
                indata["col$i"] = T[1:20;]
            end
            i = length(indata) + 1
            indata["col$i"] = [randstring(10) for j=1:20]  # ASCIIString column
            i += 1
            indata["col$i"] = ones(Bool, 20)  # Bool column
            i += 1
            indata["col$i"] = reshape([1:40;], (2, 20))  # vector Int64 column
            i += 1
            indata["col$i"] = [randstring(5) for j=1:2, k=1:20]  # vector ASCIIString col
            indata["vcol"] = [randstring(j) for j=1:20]  # variable length column
            indata["VCOL"] = [collect(1.:j) for j=1.:20.] # variable length

            # test writing
            write(f, indata; varcols=["vcol", "VCOL"])

            # test reading
            colnames = FITSIO.colnames(f[2])
            for (colname, incol) in indata
                outcol = read(f[2], colname)  # table is in extension 2 (1 = primary hdr)
                @test outcol == incol
                @test eltype(outcol) == eltype(incol)
                @test colname in colnames
            end

            # TODO: remove the tests for deprecation warnings when we don't issue them, but keep the
            # `@test_throws` test.
            @test_throws Exception @test_nowarn(read(f[2], "vcol", case_sensitive=false))
            @test_nowarn read(f[2], "col2")
            @test_logs (:warn, r"case_sensitive") read(f[2], "COL2")

            # Test representation
            @test repr(f[2])[end-38:end] == "\n\n         (*) = variable-length column"

            ## ASCII tables

            indata = Dict{String, Array}()
            for (i, T) in enumerate([Int16, Int32, Float32, Float64])
                indata["col$i"] = T[1:20;]
            end
            i = length(indata) + 1
            indata["col$i"] = [randstring(10) for j=1:20]

            write(f, indata; hdutype=ASCIITableHDU)

            # For ASCII tables, the types don't round trip so we need to define the
            # expected output type for each input type.
            expected_type = Dict(Int16=>Int32, Int32=>Int32,
                                 Float32=>Float64, Float64=>Float64,
                                 String=>String)
            colnames = FITSIO.colnames(f[3])
            for (colname, incol) in indata
                outcol = read(f[3], colname)  # table is in extension 3
                @test outcol == incol
                @test eltype(outcol) == expected_type[eltype(incol)]
                @test colname in colnames
            end
            # test show/repr on ASCIITableHDU by checking that a couple lines are what we expect
            lines = split(repr(f[3]), "\n")
            @test lines[4] == "Rows: 20"
            @test lines[6] == "         col3  Float64  E26.17  "

            # test variations on AbstractDict (issue #177)
            ordered_indata = OrderedDict(indata)
            write(f, ordered_indata; hdutype=ASCIITableHDU)

            expected_type = Dict(Int16=>Int32, Int32=>Int32,
                                 Float32=>Float64, Float64=>Float64,
                                 String=>String)
            colnames = FITSIO.colnames(f[3])
            for (colname, incol) in ordered_indata
                outcol = read(f[3], colname)  # table is in extension 3
                @test outcol == incol
                @test eltype(outcol) == expected_type[eltype(incol)]
                @test colname in colnames
            end
            # test show/repr on ASCIITableHDU by checking that a couple lines are what we expect
            lines = split(repr(f[3]), "\n")
            @test lines[4] == "Rows: 20"
            @test lines[6] == "         col3  Float64  E26.17  "
        end
    end
    @testset "Tables.jl integration" begin
        fname = tempname() * ".fits"
        FITS(fname, "w") do f
            col1 = [1., 2., 3.]
            col2 = [1, 2, 3]
            col3 = [1 2 3
                    4 5 6]
            vcol = [rand(10) for i in 1:2]
            colnames = ["col1", "col2", "col3"]
            cols = [col1, col2, col3]
            write(f, vcat(colnames, ["vcol"]), vcat(cols, [vcol]), varcols=["vcol"])
            tab = f[2]

            @testset "types work out" begin
                @test Tables.istable(FITSIO.ASCIITableHDU)
                @test Tables.istable(FITSIO.TableHDU)
                @test Tables.istable(typeof(tab))
                @test Tables.columnaccess(typeof(tab))
            end

            @test Tables.columns(tab) === tab

            @testset "columnaccess with getcolumn" begin
                @test Tables.getcolumn(tab, :col1) == col1
                @test Tables.getcolumn(tab, :col2) == col2
                @test all(Tables.getcolumn(tab, :col3) .== eachcol(col3))
                @test_throws ErrorException Tables.getcolumn(tab, :vcol)
                @test Tables.getcolumn(tab, 1) == col1
                @test Tables.getcolumn(tab, 2) == col2
                @test all(Tables.getcolumn(tab, 3) .== eachcol(col3))
                @test_throws BoundsError Tables.getcolumn(tab, 4)
            end

            @testset "column names match" begin
                @test Tables.columnnames(tab) == Symbol.(colnames)
            end

            @testset "row iteration" begin
                row = first(Tables.rows(tab))
                #@test eltype(tab) == typeof(row) #TODO should this work?
                @test row.col1 == col1[1]
                @test Tables.getcolumn(row, :col1) == col1[1]
                @test Tables.getcolumn(row, 1) == col1[1]
                @test propertynames(row) == Symbol.(colnames)
            end
        end
        rm(fname, force=true)
    end
end

@testset "FITSHeader" begin
    tempnamefits() do fname
        FITS(fname, "w") do f

            # test that show() works on an empty file and that the beginning and end
            # arre what we expect.
            s = repr(f)
            @test s[1:6] == "File: "
            @test s[end-7:end] == "No HDUs."

            @test_throws ErrorException FITSHeader(["KEY1"], [1, 2, 3], ["comment 1", "comment 2"])

            inhdr = FITSHeader(["FLTKEY", "INTKEY", "BOOLKEY", "STRKEY", "COMMENT",
                                "HISTORY"],
                               [1.0, 1, true, "string value", nothing, nothing],
                               ["floating point keyword",
                                "",
                                "boolean keyword",
                                "string value",
                                "this is a comment",
                                "this is a history"])
            inhdr2 = FITSHeader(map(NamedTuple{(:key,:value,:comment)}, [
                ("FLTKEY", 1.0, "floating point keyword",),
                ("INTKEY", 1, "",),
                ("BOOLKEY", true, "boolean keyword",),
                ("STRKEY", "string value", "string value",),
                ("COMMENT", nothing, "this is a comment",),
                ("HISTORY", nothing, "this is a history",),
            ]))
            @test inhdr.keys == inhdr2.keys
            @test inhdr.values == inhdr2.values
            @test inhdr.comments == inhdr2.comments

            @test repr(inhdr) == """
FLTKEY  =                  1.0 / floating point keyword
INTKEY  =                    1
BOOLKEY =                    T / boolean keyword
STRKEY  = 'string value'       / string value
COMMENT this is a comment
HISTORY this is a history"""

            inhdr["INTKEY"] = 2  # test setting by key
            inhdr[1] = 2.0  # test settting by index
            set_comment!(inhdr, "INTKEY", "integer keyword") # test setting a comment

            # Test reading possibly missing keyword
            @test_throws KeyError inhdr["BADKEY"]
            @test getkey(inhdr, "BADKEY", nothing) === nothing
            @test getkey(inhdr, "INTKEY", nothing) == "INTKEY"
            @test get(inhdr, "BADKEY", nothing) === nothing
            @test get(inhdr, "INTKEY", nothing) == inhdr["INTKEY"]
            @test get(() -> nothing, inhdr, "BADKEY") == nothing
            @test get(() -> nothing, inhdr, "INTKEY") == inhdr["INTKEY"]

            indata = reshape(Float32[1:100;], 5, 20)
            write(f, indata; header=inhdr)

            # Write a second block.
            inhdr2 = deepcopy(inhdr)
            inhdr2["INTKEY"] = 3 # Set it to a different value.
            write(f, indata; header=inhdr2)

            outhdr = read_header(f[1])
            @test outhdr["FLTKEY"] === 2.0
            @test outhdr["INTKEY"] === 2
            @test outhdr["BOOLKEY"] === true
            @test outhdr["STRKEY"] == "string value"
            @test get_comment(outhdr, 13) == "this is a comment"
            @test get_comment(outhdr, 14) == "this is a history"
            @test length(outhdr) == 14
            @test haskey(outhdr, "FLTKEY")

            # Read entire header as a single string
            s = read_header(f[1], String)
            @test s[1:9] == "SIMPLE  ="  # all headers should start with this.
            @test length(s) == (9 + length(inhdr)) * 80  # 9 lines = 8 default + "END"

            # Test to check that read_header gets the right block even after reading another.
            s_reread = read_header(f[1])
            s_reread = read_header(f[2])
            s_reread = read_header(f[1], String)
            @test s == s_reread

            # update an existing keyword, and read it directly
            write_key(f[1], "FLTKEY", 3.0)
            @test read_key(f[1], 9) == ("FLTKEY", 3.0, "floating point keyword")
            @test read_key(f[1], "FLTKEY") == (3.0, "floating point keyword")

            # Test appending a keyword, then modifying a keyword of different
            # values with write_key()
            for value in [1.0, "string value", 42, false, nothing]
                write_key(f[1], "NEWKEY", value, "new key comment")
                @test read_key(f[1], "NEWKEY") == (value, "new key comment")
                @test read_key(f[1], 15) == ("NEWKEY", value, "new key comment")
            end

            # Test that show() works and that the beginning of output is what we expect.
            @test repr(f)[1:6] == "File: "

            # Test the deletion of a key and verify that deleting a
            # non-existing key throws an error here.
            dhdr = deepcopy(inhdr)
            delete!(dhdr, "FLTKEY")
            @test !haskey(dhdr, "FLTKEY")

            @test_throws KeyError delete!(dhdr, "aaabbbbccccdddd")


            # Test multiple deletes
            dhdr = deepcopy(inhdr)
            delete!(dhdr, "FLTKEY")
            delete!(dhdr, "INTKEY")
            @test !haskey(dhdr, "FLTKEY") & !haskey(dhdr, "INTKEY")


        end

        hdr = FITS(fname, "r") do f
            read_header(f[1])
        end
        hdrfname = read_header(fname)
        @test keys(hdr) == keys(hdrfname)
        @test values(hdr) == values(hdrfname)
        for k in keys(hdr)
            @test get_comment(hdr, k) == get_comment(hdrfname, k)
        end
    end

    @testset "default_header" begin
        data = fill(Int16(2), 5, 6, 2)
        hdr = default_header(data)
        @test hdr isa FITSHeader
        @test hdr["SIMPLE"] == true
        @test hdr["BITPIX"] == 16
        @test hdr["NAXIS"] == 3
        @test hdr["NAXIS1"] == 2
        @test hdr["NAXIS2"] == 6
        @test hdr["NAXIS3"] == 5
        @test hdr["EXTEND"] == true
    end
end

# -----------------------------------------------------------------------------
# parsing non-standard keyword records

# `create_test_file` : Create a simple FITS file for testing, with the
# given header string added after the required keywords. The length of
# `header` must be a multiple of 80.  The purpose of creating such
# files is to test the parsing of non-standard FITS keyword records
# (non-standard files can't be created with cfitsio).

function create_test_file(fname::AbstractString, header::String)
    if length(header) % 80 != 0
        error("length of header must be multiple of 80")
    end

    open(fname, "w") do f

        stdhdr = "SIMPLE  =                    T / file does conform to FITS standard             BITPIX  =                  -64 / number of bits per data pixel                  NAXIS   =                    2 / number of data axes                            NAXIS1  =                   10 / length of data axis 1                          NAXIS2  =                   10 / length of data axis 2                          EXTEND  =                    T / FITS dataset may contain extensions            "
        endline = "END                                                                             "
        data = fill(0., (10, 10))  # 10x10 array of big-endian Float64 zeros

        # write header
        write(f, stdhdr)
        write(f, header)
        write(f, endline)

        # add padding
        block_position = (length(stdhdr) + length(header) + length(endline)) % 2880
        padding = (block_position == 0) ? 0 : 2880 - block_position
        write(f, " "^padding)

        # write data
        write(f, data)

        # add padding
        block_position = sizeof(data) % 2880
        padding = (block_position == 0) ? 0 : 2880 - block_position
        write(f, fill(0x00, (padding,)))
    end
end

@testset "Non-standard keywords" begin
    mktempdir() do dir
        fname = randomfilename(dir)
        # Create a test file with a few non-standard keyword records
        header = "        Warning: CROTA2 is inaccurate due to considerable skew                  SKEW    =  1.9511305786508E+00,  1.9234924037208E+00 /Measure of skew           "
        create_test_file(fname, header)

        # check that we can read the header (and data).
        FITS(fname) do f
            hdr = read_header(f[1])
            @test hdr["SKEW"] == "1.9511305786508E+00,"
            @test get_comment(hdr, "SKEW") == "1.9234924037208E+00 /Measure of skew"
            data = read(f[1])
            @test data == zeros(10, 10)
        end

        # Test that we can read is at a memory backed file
        FITS(read(fname)) do f
            hdr = read_header(f[1])
            data = read(f[1])
            @test data == zeros(10, 10)
            # test that we can write it back out.
            FITS(randomfilename(dir), "w") do f2
                write(f2, data; header=hdr)
            end
        end
    end
end
@testset "Variable Length Column Operations" begin
    # Helper function to create test data
    function create_test_data()
        return (
            strings = ["short", "medium length", "very long string"],
            numbers = [[1.0, 2.0], [3.0, 4.0, 5.0], [6.0]],
            integers = [[1], [1,2,3], [1,2]]
        )
    end

    tempnamefits() do fname
        # Test writing variable length columns
        let
            test_data = create_test_data()
            
            FITS(fname, "w") do f
                # Write primary HDU first
                write(f, zeros(Float32, 1))

                # Create data dictionary
                data = Dict(
                    "STRINGS" => test_data.strings,
                    "NUMBERS" => test_data.numbers,
                    "INTEGERS" => test_data.integers
                )

                # Write table with variable length columns
                write(f, data; varcols=["STRINGS", "NUMBERS", "INTEGERS"])
            end
        end

        # Test reading variable length columns
        let
            test_data = create_test_data()
            
            FITS(fname, "r") do f
                hdu = f[2]  # Skip primary HDU

                # Test string column
                @test read(hdu, "STRINGS") == test_data.strings

                # Test numeric columns with explicit type conversion
                numbers = read(hdu, "NUMBERS")
                @test length(numbers) == length(test_data.numbers)
                for (i, (orig, read_val)) in enumerate(zip(test_data.numbers, numbers))
                    @test length(read_val) == length(orig)
                    @test all(read_val .== orig)
                end

                integers = read(hdu, "INTEGERS")
                @test length(integers) == length(test_data.integers)
                for (i, (orig, read_val)) in enumerate(zip(test_data.integers, integers))
                    @test length(read_val) == length(orig)
                    @test all(read_val .== orig)
                end
            end
        end
    end
end

@testset "Error Handling for Variable Length Columns" begin
    tempnamefits() do fname
        FITS(fname, "w") do f
            # Write primary HDU
            write(f, zeros(Float32, 1))

            # Test ASCII table restriction
            @test_throws ErrorException begin
                write(f, Dict("TEST" => [[1,2], [3,4]]); 
                      hdutype=ASCIITableHDU, 
                      varcols=["TEST"])
            end

            # Test unsupported type
            @test_throws ErrorException begin
                complex_data = [[1+2im], [3+4im]]
                write(f, Dict("COMPLEX" => complex_data); 
                      varcols=["COMPLEX"])
            end
        end
    end
end

@testset "Table Writing" begin
    tempnamefits() do fname
        FITS(fname, "w") do f
            # Write primary HDU
            write(f, zeros(Float32, 1))

            # Create simple table
            data = Dict(
                "COL1" => [1, 2, 3],
                "COL2" => ["a", "b", "c"]
            )

            # Write table
            write(f, data)

            # Verify data
            hdu = f[2]
            @test read(hdu, "COL1") == [1, 2, 3]
            @test read(hdu, "COL2") == ["a", "b", "c"]
        end
    end
end

@testset "Variable Column Type Handling" begin
    tempnamefits() do fname
        FITS(fname, "w") do f
            # Write primary HDU
            write(f, zeros(Float32, 1))

            # Test different variable length types
            data = Dict(
                "VAR_INT" => [[1], [1,2], [1,2,3]],
                "VAR_FLOAT" => [[1.0], [1.0,2.0], [1.0,2.0,3.0]],
                "VAR_STRING" => ["a", "bb", "ccc"]
            )

            # Write data
            write(f, data; varcols=["VAR_INT", "VAR_FLOAT", "VAR_STRING"])

            # Verify data
            hdu = f[2]
            @test all(read(hdu, "VAR_INT")[i] == data["VAR_INT"][i] 
                     for i in 1:3)
            @test all(read(hdu, "VAR_FLOAT")[i] == data["VAR_FLOAT"][i] 
                     for i in 1:3)
            @test read(hdu, "VAR_STRING") == data["VAR_STRING"]
        end
    end
    @testset "Error Handling for Invalid Types" begin
        tempnamefits() do fname
            FITS(fname, "w") do f
                GC.@preserve begin
                    write(f, zeros(Float32, 1))
                    
                    # Test deeply nested vector rejection
                    nested_data = Dict("NESTED" => [[[1,2], [3,4]]])
                    @test_throws ErrorException write(f, nested_data; varcols=["NESTED"])
                end
            end
        end
    end
end
# Define test data creation function
function create_test_data()
    strings = ["short", "medium length", "very long string"]
    numbers = Vector{Float64}[[1.0, 2.0], [3.0, 4.0, 5.0], [6.0]]
    integers = Vector{Int64}[[1], [1,2,3], [1,2]]
    return strings, numbers, integers
end

@testset "Mixed Variable Length Types" begin
    tempnamefits() do fname
        strings, numbers, integers = create_test_data()

        FITS(fname, "w") do f
            # Ensure primary HDU is created properly
            GC.@preserve f strings numbers integers begin
                write(f, zeros(Float32, 1))  # Primary HDU
                
                # Test writing multiple variable length columns of different types together
                data = Dict(
                    "STRINGS" => strings,
                    "NUMBERS" => numbers,
                    "INTEGERS" => integers
                )
                write(f, data; varcols=["STRINGS", "NUMBERS", "INTEGERS"])
            end
        end

        FITS(fname, "r") do f
            hdu = f[2]
            
            # Test reading and verifying mixed types with GC protection
            GC.@preserve f hdu strings numbers integers begin
                read_strings = read(hdu, "STRINGS")
                read_numbers = read(hdu, "NUMBERS")
                read_integers = read(hdu, "INTEGERS")
                
                # Version-specific string comparison
                @static if VERSION < v"1.6"
                    # Direct string comparison for older Julia
                    @test all(String.(read_strings) .== String.(strings))
                else
                    # Modern comparison
                    @test read_strings == strings
                end
                
                # Verify numeric data
                @test all(length.(read_numbers) .== length.(numbers))
                @test all(length.(read_integers) .== length.(integers))
                
                # Verify content equality with GC protection
                for (orig, read) in zip(numbers, read_numbers)
                    GC.@preserve orig read begin
                        @test all(orig .== read)
                    end
                end
                
                for (orig, read) in zip(integers, read_integers)
                    GC.@preserve orig read begin
                        @test all(orig .== read)
                    end
                end
            end
        end
    end
end

if Sys.iswindows()
    ENV["JULIA_CPU_THREADS"] = "1"
    Base.GC.enable(true)
end

@testset "Variable Length Column Operations" begin
    @testset "Reading Operations" begin
        tempnamefits() do fname
            FITS(fname, "w") do f
                # Write primary HDU first
                write(f, zeros(Float32, 1))

                # Prepare test data
                let
                    string_data = ["short", "medium length", "very long string"]
                    numeric_data = [[1.0, 2.0], [3.0, 4.0, 5.0], [6.0]]
                    int_data = [[1], [1,2,3], [1,2]]

                    # Write data
                    data = Dict(
                        "STRINGS" => string_data,
                        "NUMBERS" => numeric_data,
                        "INTEGERS" => int_data
                    )
                    write(f, data; varcols=["STRINGS", "NUMBERS", "INTEGERS"])
                end
            end

            # Read and verify in separate FITS instance
            FITS(fname, "r") do f
                hdu = f[2]  # Skip primary HDU
                
                # Test string column
                let
                    expected = ["short", "medium length", "very long string"]
                    result = read(hdu, "STRINGS")
                    @test result == expected
                end

                # Test numeric column
                let
                    expected = [[1.0, 2.0], [3.0, 4.0, 5.0], [6.0]]
                    result = read(hdu, "NUMBERS")
                    @test length(result) == length(expected)
                    for (r, e) in zip(result, expected)
                        @test r == e
                    end
                end

                # Test integer column
                let
                    expected = [[1], [1,2,3], [1,2]]
                    result = read(hdu, "INTEGERS")
                    @test length(result) == length(expected)
                    for (r, e) in zip(result, expected)
                        @test r == e
                    end
                end
            end
        end
    end

    @testset "Writing Operations" begin
        tempnamefits() do fname
            FITS(fname, "w") do f
                # Write primary HDU
                write(f, zeros(Float32, 1))

                # Test string column
                let
                    data = ["short", "medium length", "very long string"]
                    write(f, Dict("STRINGS" => data); varcols=["STRINGS"])
                    result = read(f[end], "STRINGS")
                    @test result == data
                end

                # Test numeric column
                let
                    data = [[1.0, 2.0], [3.0, 4.0, 5.0], [6.0]]
                    write(f, Dict("NUMBERS" => data); varcols=["NUMBERS"])
                    result = read(f[end], "NUMBERS")
                    @test all(result[i] == data[i] for i in 1:length(data))
                end

                # Test error cases
                @test_throws ErrorException write(f, Dict("TEST" => [1,2,3]); 
                    hdutype=ASCIITableHDU, 
                    varcols=["TEST"])

                @test_throws ErrorException write(f, 
                    Dict("COMPLEX" => [[1+2im], [3+4im]]); 
                    varcols=["COMPLEX"])
            end
        end
    end

    @testset "Table Functions" begin
        tempnamefits() do fname
            FITS(fname, "w") do f
                # Write primary HDU
                write(f, zeros(Float32, 1))

                let
                    # Basic table
                    colnames = ["COL1", "COL2"]
                    coldata = [[1, 2, 3], ["a", "b", "c"]]
                    write(f, Dict(zip(colnames, coldata)))
                    
                    @test read(f[end], "COL1") == [1, 2, 3]
                    @test read(f[end], "COL2") == ["a", "b", "c"]

                    # Test with units
                    units = Dict("COL1" => "m/s", "COL2" => "kg")
                    write(f, Dict(zip(colnames, coldata)); units=units)

                    # Test with version
                    write(f, Dict(zip(colnames, coldata)); ver=2)
                    @test read_key(f[end], "EXTVER") == (2, "")

                    # Test with header
                    hdr = FITSHeader(["TEST"], [1], ["test comment"])
                    write(f, Dict(zip(colnames, coldata)); header=hdr)
                    @test read_key(f[end], "TEST") == (1, "test comment")

                    # Test variable length
                    var_data = [[[1], [1,2]], ["a", "bb"]]
                    write(f, Dict(zip(colnames, var_data)); varcols=["COL1"])
                    result = read(f[end], "COL1")
                    @test all(result[i] == var_data[1][i] for i in 1:length(var_data[1]))
                end
            end
        end
    end
end
@test_deprecated FITSIO.libcfitsio_version()

# test we can still access Libcfitsio

@test begin
    using FITSIO.Libcfitsio
    Libcfitsio.libcfitsio_version() isa VersionNumber
end
