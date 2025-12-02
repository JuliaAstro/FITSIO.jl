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
            @test get(() -> nothing, inhdr, "BADKEY") === nothing
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
