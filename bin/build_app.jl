# Run with Julia 1.13.2 and PackageCompiler available in the active build environment:
# julia --project=/path/to/build-env bin/build_app.jl /path/to/new-bundle
using Pkg
using PackageCompiler

VERSION==v"1.13.2" || error("Build myloc with Julia 1.13.2; current Julia is $VERSION")

length(ARGS)==1 || error("Usage: julia --project=BUILD_ENV bin/build_app.jl NEW_BUNDLE_DIRECTORY")
root=dirname(@__DIR__)
destination=abspath(only(ARGS))
ispath(destination) && error("Bundle destination already exists: $destination")
mktempdir() do staging
    # Only runtime source enters the compilation project. Build dependencies
    # and local research files are not runtime dependencies or release content.
    cp(joinpath(root,"Project.toml"),joinpath(staging,"Project.toml"))
    cp(joinpath(root,"src"),joinpath(staging,"src"))
    cp(joinpath(root,"ext"),joinpath(staging,"ext"))
    Pkg.activate(staging)
    Pkg.instantiate()
    # Exercise build-only calculations while collecting code;
    # PackageCompiler still builds the image itself on its required one thread.
    withenv("JULIA_NUM_THREADS"=>"2") do
        create_app(staging,destination;executables=["myloc"=>"julia_main"],
            precompile_execution_file=joinpath(@__DIR__,"precompile_app.jl"),
            cpu_target="generic",incremental=false,filter_stdlibs=false)
    end
    cp(normpath(Sys.BINDIR,"..","LICENSE.md"),joinpath(destination,"Julia-LICENSE.md"))
    if Sys.iswindows()
        powershell=joinpath(ENV["SystemRoot"],"System32","WindowsPowerShell","v1.0","powershell.exe")
        run(`$powershell -NoProfile -File $(joinpath(@__DIR__,"add_manifest.ps1")) -Executable $(joinpath(destination,"bin","myloc.exe"))`)
    end
    share=joinpath(destination,"share","maoyeloc")
    mkpath(share)
    cp(joinpath(root,"src","data"),joinpath(share,"data"))
    # Keep the public numerical source referenced by the mathematical manual.
    cp(joinpath(root,"src"),joinpath(destination,"src"))
    cp(joinpath(root,"ext"),joinpath(destination,"ext"))
    cp(joinpath(root,"Project.toml"),joinpath(destination,"Project.toml"))
    for name in ("maoyeloc.jl","myloc.jl","batch.jl","plot_locations.jl","review_catalog.jl","generate_travel_time_table.jl")
        cp(joinpath(@__DIR__,name),joinpath(destination,"bin",name))
    end
    cp(joinpath(root,"README.md"),joinpath(destination,"README.md"))
    mkpath(joinpath(destination,"docs"))
    for name in ("用户手册.md","数据约定.md","MaoyeLoc原理.md")
        cp(joinpath(root,"docs",name),joinpath(destination,"docs",name))
    end
    if Sys.iswindows()
        cp(joinpath(@__DIR__,"install.ps1"),joinpath(destination,"install.ps1"))
        cp(joinpath(@__DIR__,"myloc.cmd"),joinpath(destination,"myloc.cmd"))
        write(joinpath(destination,"install.cmd"),"@echo off\r\npowershell -NoProfile -ExecutionPolicy Bypass -File \"%~dp0install.ps1\"\r\npause\r\n")
    else
        cp(joinpath(@__DIR__,"install.sh"),joinpath(destination,"install.sh"))
    end
end
println("Built ",joinpath(destination,"bin",Sys.iswindows() ? "myloc.exe" : "myloc"))
