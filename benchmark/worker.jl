# One benchmark point (case, size) in a fresh process; called by run.jl.
# usage: julia --project=<env> worker.jl <case> <size> <csvfile>

using FewBodyToolkit, LinearAlgebra, Printf
include("cases.jl")

# peak memory: VmHWM of this process, reset by writing 5 to clear_refs (Linux)
procstat(key) = parse(Int, split(only(filter(l -> startswith(l, key), readlines("/proc/self/status"))))[2]) * 1024
function measure(f)
    GC.gc(); GC.gc()
    write("/proc/self/clear_refs", "5")
    rss0 = procstat("VmRSS:")
    t = @elapsed y = f()
    return y, t, procstat("VmHWM:") - rss0
end

const COLUMNS = ["case","n","N","t_solve","t_build","t_inv1","t_inv6","t_inv6_prep",
                 "dpeak_solve","dpeak_build","dpeak_inv1","rss0","E1_re","E1_im","E2_re","E2_im","E3_re","E3_im"]

function run_point(case, n)
    r = Dict{String,Any}(c => NaN for c in COLUMNS)
    r["case"] = case.name; r["n"] = n
    pp, np = params(case.kind, n)
    solve, matrices = (isdefined(FewBodyToolkit, f) ? getfield(FewBodyToolkit, f) : nothing for f in solvers(case.kind))
    if case.kind == :inverse
        (T, V, S), r["t_build"], r["dpeak_build"] = measure(() -> ISGL_matrices(pp, np))
        lam, r["t_inv1"], r["dpeak_inv1"] = measure(() -> inverse_solve(T, V, S; target_energy=0.0, threshold=np.threshold))
        r["t_inv6"] = @elapsed foreach(e -> inverse_solve(T, V, S; target_energy=e, threshold=np.threshold), INVERSE_ENERGIES)
        if isdefined(FewBodyToolkit, :inverse_prepare)
            r["t_inv6_prep"] = @elapsed begin
                P = FewBodyToolkit.inverse_prepare(T, V, S; threshold=np.threshold)
                foreach(e -> inverse_solve(P; target_energy=e), INVERSE_ENERGIES)
            end
        end
        r["N"] = size(S, 1)
        E = lam
    else
        E, r["t_solve"], r["dpeak_solve"] = measure(() -> solve(pp, np; solve_kw(case.kw)...))
        r["N"] = length(E) # one entry per basis function
        if !isnothing(matrices) # >= 0.5.4
            M, r["t_build"], r["dpeak_build"] = measure(() -> matrices(pp, np; case.kw...))
            M = nothing
        end
    end
    E = sort(E[abs.(E) .> 0], by=real) # entries cut by the overlap threshold are 0
    for i in 1:3
        r["E$(i)_re"] = real(E[i]); r["E$(i)_im"] = imag(E[i])
    end
    GC.gc(); r["rss0"] = procstat("VmRSS:")
    return r
end

case = only(filter(c -> c.name == ARGS[1], CASES))
run_point(case, WARMUP_SIZE) # compilation
r = run_point(case, parse(Int, ARGS[2]))
open(ARGS[3], "a") do io
    filesize(ARGS[3]) == 0 && println(io, join(COLUMNS, ","))
    println(io, join([r[c] isa AbstractFloat ? @sprintf("%.10g", r[c]) : string(r[c]) for c in COLUMNS], ","))
end
@printf("%-11s n=%3d N=%6d  t_solve=%8.2f s  t_build=%8.2f s  dpeak_solve=%7.3f GB\n",
        r["case"], r["n"], r["N"], r["t_solve"], r["t_build"], r["dpeak_solve"]/1e9)
