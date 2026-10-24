# Runs all benchmark cases (cases.jl) for one FewBodyToolkit version, each point in a fresh process.
# usage: julia [-t <n>] benchmark/run.jl <fbtk> [--quick] [--cases=name1,name2]
#   <fbtk>: a registered version (e.g. 0.5.4) or the path to a FewBodyToolkit.jl checkout
#   --quick: only the smallest size of each case (to check the setup)
# Output in benchmark/results/: <label>.csv (one row per point), <label>.toml (metadata),
# <label>_Manifest.toml (environment), <label>.log. <label> = <version or commit>_<date>.
# No host name or local paths are recorded (results may be published).

using Pkg, TOML, Dates, LinearAlgebra

const HERE = @__DIR__
const OPTS = filter(a -> startswith(a, "--"), ARGS)
const SPEC = only(filter(a -> !startswith(a, "--"), ARGS))
const QUICK = "--quick" in OPTS
const CASEFILTER = [split(o[9:end], ",") for o in OPTS if startswith(o, "--cases=")]

include("cases.jl") # only definitions; the package itself is loaded in the workers

# environment with the requested version
isrepo = isdir(SPEC)
env = mktempdir(; cleanup=false)
Pkg.activate(env; io=devnull)
isrepo ? Pkg.develop(path=abspath(SPEC); io=devnull) : Pkg.add(name="FewBodyToolkit", version=SPEC; io=devnull)
Pkg.precompile()
manifest = TOML.parsefile(joinpath(env, "Manifest.toml"))
fbtk = only(manifest["deps"]["FewBodyToolkit"])
gitcmd(args...) = readchomp(`git -C $(abspath(SPEC)) $args`)
if isrepo
    commit = gitcmd("rev-parse", "--short", "HEAD")
    dirty = !isempty(gitcmd("status", "--porcelain", "--untracked-files=no"))
    id = "commit-" * commit * (dirty ? "-dirty" : "")
else
    id = "v" * fbtk["version"]
end
label = id * "_" * Dates.format(now(), "yyyymmdd-HHMM") * (QUICK ? "_quick" : "")
resdir = joinpath(HERE, "results"); mkpath(resdir)
csvfile = joinpath(resdir, label * ".csv")
write(joinpath(resdir, label * "_Manifest.toml"),  # without the local path of a developed checkout
      replace(read(joinpath(env, "Manifest.toml"), String), r"\npath = \".*\"" => "\npath = \"<local checkout>\""))

meta = Dict(
    "label" => label, "date" => string(now()),
    "fbtk_spec" => isrepo ? "local checkout" : SPEC, "fbtk_version" => fbtk["version"],
    "fbtk_commit" => isrepo ? gitcmd("rev-parse", "HEAD") : get(fbtk, "git-tree-sha1", ""),
    "fbtk_dirty" => isrepo ? dirty : false,
    "julia_version" => string(VERSION), "cpu" => Sys.cpu_info()[1].model, "logical_cores" => Sys.CPU_THREADS,
    "julia_threads" => Threads.nthreads(), "blas_threads" => BLAS.get_num_threads(), "blas" => string(BLAS.get_config()),
    "loadavg_start" => read("/proc/loadavg", String) |> chomp,
    "quick" => QUICK,
)
open(io -> TOML.print(io, meta), joinpath(resdir, label * ".toml"), "w")

worker = joinpath(HERE, "worker.jl")
open(joinpath(resdir, label * ".log"), "w") do log
    for case in CASES
        isempty(CASEFILTER) || case.name in only(CASEFILTER) || continue
        VersionNumber(fbtk["version"]) >= case.since || continue
        for n in (QUICK ? case.sizes[1:1] : case.sizes)
            cmd = `$(Base.julia_cmd()) --project=$env --threads=$(Threads.nthreads()) $worker $(case.name) $n $csvfile`
            println(log, now(), "  ", case.name, " n=", n, "  load ", chomp(read("/proc/loadavg", String))); flush(log)
            ok = success(pipeline(cmd; stdout=log, stderr=log))
            ok || println(log, "FAILED (an empty or missing row with no error message usually means out of memory)")
            flush(log)
        end
    end
    println(log, now(), "  done")
end
println("results: ", csvfile)
