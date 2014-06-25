using BinDeps
@BinDeps.setup

version = "3360"
url = "ftp://heasarc.gsfc.nasa.gov/software/fitsio/c/cfitsio$version.tar.gz"

libcfitsio = library_dependency("libcfitsio")
provides(Sources, URI(url), libcfitsio, unpacked_dir="cfitsio")
depsdir = BinDeps.depsdir(libcfitsio)
srcdir = joinpath(depsdir, "src", "cfitsio")
prefix = joinpath(depsdir, "usr")
@unix_only libfilename = "libcfitsio.so"
@osx_only libfilename = "libcfitsio.dylib"
provides(BuildProcess,
         (@build_steps begin
            GetSources(libcfitsio)
            @build_steps begin
                ChangeDirectory(srcdir)
                FileRule(joinpath(prefix,"lib",libfilename),
                         @build_steps begin
                         `./configure --prefix=$prefix`
                         `make shared install`
                         end)
            end
         end),
         libcfitsio)

@BinDeps.install [:libcfitsio => :libcfitsio]
