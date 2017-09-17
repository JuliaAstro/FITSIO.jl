using BinDeps
using Compat
using Compat.Sys: iswindows, isapple, isunix

@BinDeps.setup

version = "3370"
baseurl = "ftp://heasarc.gsfc.nasa.gov/software/fitsio/c/"

if isunix()
    archivename = "cfitsio$(version).tar.gz"
elseif iswindows()
    archivename = "cfitsio_MSVC_$(Sys.WORD_SIZE)bit_DLL_$(version).zip"
end

libcfitsio = library_dependency("libcfitsio", aliases=["cfitsio"])
downloadsdir = BinDeps.downloadsdir(libcfitsio)
libdir = BinDeps.libdir(libcfitsio)
srcdir = BinDeps.srcdir(libcfitsio)

if isapple()
    libfilename = "libcfitsio.dylib"
elseif isunix()
    libfilename = "libcfitsio.so"
elseif iswindows()
    libfilename = "cfitsio.dll"
end

# Unix
prefix = joinpath(BinDeps.depsdir(libcfitsio), "usr")
provides(Sources, URI(baseurl*archivename), libcfitsio, unpacked_dir="cfitsio")
provides(BuildProcess,
         (@build_steps begin
             GetSources(libcfitsio)
             @build_steps begin
                 ChangeDirectory(joinpath(srcdir, "cfitsio"))
                 FileRule(joinpath(libdir, libfilename),
                          @build_steps begin
                              `./configure --prefix=$prefix --enable-reentrant`
                              `make shared install`
                          end)
             end
          end), libcfitsio, os = :Unix)

# Windows
provides(BuildProcess,
	 (@build_steps begin
             FileDownloader(baseurl*archivename,
                            joinpath(downloadsdir, archivename))
	     CreateDirectory(srcdir, true)
	     FileUnpacker(joinpath(downloadsdir, archivename), srcdir,
                          joinpath(srcdir, libfilename))
	     CreateDirectory(libdir, true)
	     @build_steps begin
		 ChangeDirectory(srcdir)
	         FileRule(joinpath(libdir, libfilename), @build_steps begin
		          `powershell -Command "cp $(libfilename) $(joinpath(libdir, libfilename))"`
			  end)
             end
	  end), libcfitsio, os = :Windows)

if iswindows()
    push!(BinDeps.defaults, BuildProcess)
end

# OSX
if isapple()
    using Homebrew
    provides(Homebrew.HB, "cfitsio", libcfitsio, os=:Darwin)
end

@BinDeps.install Dict(:libcfitsio => :libcfitsio)

if iswindows()
    pop!(BinDeps.defaults)
end
