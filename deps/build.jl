using BinDeps
using Compat

@BinDeps.setup

version = "3370"
baseurl = "ftp://heasarc.gsfc.nasa.gov/software/fitsio/c/"

if is_unix()
    archivename = "cfitsio$(version).tar.gz"
elseif is_windows()
    archivename = "cfitsio_MSVC_$(Sys.WORD_SIZE)bit_DLL_$(version).zip"
end

libcfitsio = library_dependency("libcfitsio", aliases=["cfitsio"])
downloadsdir = BinDeps.downloadsdir(libcfitsio)
libdir = BinDeps.libdir(libcfitsio)
srcdir = BinDeps.srcdir(libcfitsio)

if is_apple()
    libfilename = "libcfitsio.dylib"
elseif is_unix()
    libfilename = "libcfitsio.so"
elseif is_windows()
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

if is_windows()
    push!(BinDeps.defaults, BuildProcess)
end

@BinDeps.install @compat Dict(:libcfitsio => :libcfitsio)

if is_windows()
    pop!(BinDeps.defaults)
end
