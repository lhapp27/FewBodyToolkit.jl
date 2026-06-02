# ISGL-only hotspot validation benchmark

using FewBodyToolkit
using BenchmarkTools
using Printf

const IS = FewBodyToolkit.ISGL
const COMPLEX_SCALING = false
const RETURN_WAVEFUNCTIONS = false
const OBS_PARAMS = (
    stateindices = Int[],
    centobs_arr = [Any[] for _ in 1:3],
    R2_arr = Int[0, 0, 0],
)

function make_case(; nmax::Int, kmax_interpol::Int, lmax::Int = 0, Lmax::Int = 0)
    vg(r) = -10.0 * exp(-r^2)
    pp = FewBodyToolkit.make_phys_params3B3D(
        masses = [1.0, 2.0, 3.0],
        species = [:x, :y, :z],
        interactions = [[vg], [vg], [vg]],
    )

    gp = (; nmax = nmax, Nmax = nmax, r1 = 0.1, rnmax = 100.0, R1 = 0.1, RNmax = 100.0)
    np = FewBodyToolkit.make_num_params3B3D(
        gem_params = gp,
        kmax_interpol = kmax_interpol,
        lmin = 0,
        lmax = lmax,
        Lmin = 0,
        Lmax = Lmax,
        threshold = 1e-10,
    )

    return pp, np
end

function median_stats(trial)
    m = median(trial)
    return (time_ms = m.time / 1.0e6, memory_mib = m.memory / 1024^2, allocs = m.allocs)
end

function print_table(title::String, rows)
    println("\n" * title)
    println(repeat("-", length(title)))
    @printf("%-22s %14s %14s %14s\n", "Stage", "Time [ms]", "Memory [MiB]", "Allocs")
    @printf("%-22s %14s %14s %14s\n", repeat("-", 22), repeat("-", 14), repeat("-", 14), repeat("-", 14))

    total_ms = 0.0
    for row in rows
        total_ms += row.time_ms
        @printf("%-22s %14.3f %14.3f %14d\n", row.name, row.time_ms, row.memory_mib, row.allocs)
    end

    println()
    @printf("Total isolated stage time: %.3f ms\n", total_ms)
end

function benchmark_stages(pp, np; samples::Int = 8)
    obs = OBS_PARAMS
    hbar_val = pp.hbar

    trial_sanity = @benchmark IS.sanity_checks($pp) samples = samples evals = 1

    trial_size = @benchmark IS.size_estimate($pp, $np, $obs, $COMPLEX_SCALING) samples = samples evals = 1

    trial_prealloc = @benchmark IS.preallocate_data($pp, $np, $obs, size_params, $COMPLEX_SCALING) setup = begin
        size_params = IS.size_estimate($pp, $np, $obs, $COMPLEX_SCALING)
    end samples = samples evals = 1

    trial_precompute = @benchmark IS.precompute_ISGL($pp, $np, size_params, precomp_arrs, temp_arrs) setup = begin
        size_params = IS.size_estimate($pp, $np, $obs, $COMPLEX_SCALING)
        precomp_arrs, temp_arrs, interpol_arrs, fill_arrs, result_arrs = IS.preallocate_data($pp, $np, $obs, size_params, $COMPLEX_SCALING)
    end samples = samples evals = 1

    trial_interpol = @benchmark IS.interpolNshoulder($pp, $np, $obs, size_params, precomp_arrs, interpol_arrs, $RETURN_WAVEFUNCTIONS, $COMPLEX_SCALING) setup = begin
        size_params = IS.size_estimate($pp, $np, $obs, $COMPLEX_SCALING)
        precomp_arrs, temp_arrs, interpol_arrs, fill_arrs, result_arrs = IS.preallocate_data($pp, $np, $obs, size_params, $COMPLEX_SCALING)
        IS.precompute_ISGL($pp, $np, size_params, precomp_arrs, temp_arrs)
    end samples = samples evals = 1

    trial_fill = @benchmark IS.fill_TVS($np, size_params, precomp_arrs, interpol_arrs, fill_arrs, $COMPLEX_SCALING, $hbar_val, false) setup = begin
        size_params = IS.size_estimate($pp, $np, $obs, $COMPLEX_SCALING)
        precomp_arrs, temp_arrs, interpol_arrs, fill_arrs, result_arrs = IS.preallocate_data($pp, $np, $obs, size_params, $COMPLEX_SCALING)
        IS.precompute_ISGL($pp, $np, size_params, precomp_arrs, temp_arrs)
        IS.interpolNshoulder($pp, $np, $obs, size_params, precomp_arrs, interpol_arrs, $RETURN_WAVEFUNCTIONS, $COMPLEX_SCALING)
    end samples = samples evals = 1

    trial_solve = @benchmark IS.solveHS($np, fill_arrs, result_arrs, $RETURN_WAVEFUNCTIONS) setup = begin
        size_params = IS.size_estimate($pp, $np, $obs, $COMPLEX_SCALING)
        precomp_arrs, temp_arrs, interpol_arrs, fill_arrs, result_arrs = IS.preallocate_data($pp, $np, $obs, size_params, $COMPLEX_SCALING)
        IS.precompute_ISGL($pp, $np, size_params, precomp_arrs, temp_arrs)
        IS.interpolNshoulder($pp, $np, $obs, size_params, precomp_arrs, interpol_arrs, $RETURN_WAVEFUNCTIONS, $COMPLEX_SCALING)
        IS.fill_TVS($np, size_params, precomp_arrs, interpol_arrs, fill_arrs, $COMPLEX_SCALING, $hbar_val, false)
    end samples = samples evals = 1

    names = [
        "sanity_checks",
        "size_estimate",
        "preallocate_data",
        "precompute_ISGL",
        "interpolNshoulder",
        "fill_TVS",
        "solveHS",
    ]
    trials = [trial_sanity, trial_size, trial_prealloc, trial_precompute, trial_interpol, trial_fill, trial_solve]

    rows = NamedTuple[]
    for (name, trial) in zip(names, trials)
        s = median_stats(trial)
        push!(rows, (name = name, time_ms = s.time_ms, memory_mib = s.memory_mib, allocs = s.allocs))
    end

    return rows
end

function parse_nmaxs(arg::String)
    parts = split(arg, ",")
    vals = Int[]
    for p in parts
        push!(vals, parse(Int, strip(p)))
    end
    return vals
end

function parse_int_list(arg::String)
    parts = split(arg, ",")
    vals = Int[]
    for p in parts
        push!(vals, parse(Int, strip(p)))
    end
    return vals
end

function run_hotspot(; nmaxs::Vector{Int} = [10, 16], samples::Int = 2, kmax_interpol::Int = 1000, lmaxs::Vector{Int} = [0], Lmax::Int = 0)
    println("ISGL hotspot validation benchmark")
    @printf("nmax values: %s | lmax values: %s | samples: %d | kmax_interpol: %d | Lmax: %d\n", string(nmaxs), string(lmaxs), samples, kmax_interpol, Lmax)

    for lmax in lmaxs
        for nmax in nmaxs
            pp, np = make_case(nmax = nmax, kmax_interpol = kmax_interpol, lmax = lmax, Lmax = Lmax)
            rows = benchmark_stages(pp, np; samples = samples)
            print_table("nmax = $(nmax), lmax = $(lmax)", rows)
        end
    end
end

function main()
    nmaxs = isempty(ARGS) ? [10, 16] : parse_nmaxs(ARGS[1])
    samples = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 2
    kmax_interpol = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 1000
    lmaxs = length(ARGS) >= 4 ? parse_int_list(ARGS[4]) : [0]
    Lmax = length(ARGS) >= 5 ? parse(Int, ARGS[5]) : 0

    run_hotspot(; nmaxs = nmaxs, samples = samples, kmax_interpol = kmax_interpol, lmaxs = lmaxs, Lmax = Lmax)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
