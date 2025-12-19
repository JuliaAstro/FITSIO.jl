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

            @test_throws BoundsError f[100]

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
        a = ones(3,3)
        tempnamefits() do fname
            FITS(fname, "w") do f
                write(f, a, name="a")
                @test read(f["a"]) == a
                @test read(f[1])   == a
                @test haskey(f, "a")
                @test !haskey(f, "b")
                @test haskey(f, 1)
                @test !haskey(f, 2)
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
