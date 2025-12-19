import Tables

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
            codes = Dict("Int16" => "I7", "Int32" => "I12",
                         "Float32" => "E26.17", "Float64" => "E26.17",
                         "String" => "A10")
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
            lines = strip.(lines[6:end])
            # we sort the lines as the order of the columns is not guaranteed
            sort!(lines, by = x -> split(x, " ")[1])
            lines_tokens = split.(lines, " ", keepempty=false)
            function tokenvec(k, d)
                el = eltype(d[k])
                el_exp = expected_type[el]
                [k, string(el_exp), codes[string(el)]]
            end
            for i in eachindex(lines_tokens)
                @test lines_tokens[i] == tokenvec("col$i", indata)
            end

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
            lines = strip.(lines[6:end])
            lines_tokens = split.(lines, " ", keepempty=false)

            keys_ord_vec = collect(keys(ordered_indata))
            # check that the order is preserved
            for (i, k) in enumerate(keys_ord_vec)
                @test lines_tokens[i] == tokenvec(k, ordered_indata)
            end
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

    @testset "ensure correct HDU pointer in indexing" begin
        tempnamefits() do fname
            FITS(fname, "w") do f
                write(f, Dict("col1" => [1, 2, 3], "col2" => [4, 5, 6]))
                @test f[1] isa FITSIO.ImageHDU{<:Any,0}
                @test f[2] isa FITSIO.TableHDU
            end
        end
    end
end
