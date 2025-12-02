using ParallelTestRunner: runtests, find_tests, parse_args
using FITSIO

const init_code = quote
    using FITSIO
    using CFITSIO
    using OrderedCollections

    # Deal with compatibility issues.
    using Test
    using Random # for `randstring`

    randomfilename(dir) = joinpath(dir, randstring(12)*".fits")
    function tempnamefits(f)
        mktempdir() do dir
            f(randomfilename(dir))
        end
    end
end

args = parse_args(Base.ARGS)
testsuite = find_tests(@__DIR__)

runtests(FITSIO, args; testsuite, init_code)
