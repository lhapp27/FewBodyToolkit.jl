# Compares benchmark result files (run.jl): times, peak memory, and the lowest energies.
# usage: julia benchmark/compare.jl results/<old>.csv results/<new>.csv     (detailed comparison)
#        julia benchmark/compare.jl results/<a>.csv results/<b>.csv ...    (overview of several versions)

using DelimitedFiles, Printf

function load(file)
    d, h = readdlm(file, ','; header=true)
    rows = Dict{Tuple{String,Int},Dict{String,Any}}()
    for i in axes(d, 1)
        rows[(string(d[i,1]), Int(d[i,2]))] = Dict(string(h[j]) => d[i,j] for j in axes(d, 2))
    end
    return rows
end

ratio(a, b) = isfinite(a) && isfinite(b) && a > 0 ? @sprintf("%6.2f", b/a) : "     -"
gb(x) = isfinite(x) ? @sprintf("%7.3f", x/1e9) : "      -"
sec(x) = isfinite(x) ? @sprintf("%8.2f", x) : "       -"

# c1 of dpeak = c0 + c1*N^2 (least squares), in bytes per N^2
function fit_c1(N, y)
    ok = isfinite.(y)
    count(ok) >= 2 || return NaN
    A = hcat(ones(count(ok)), float.(N[ok]).^2)
    return (A \ y[ok])[2]
end

# overview: time and peak memory of the solve (inverse: one inverse_solve) per version, and the
# largest deviation of the lowest energies from the last file
if length(ARGS) > 2
    R = load.(ARGS)
    names = [first(split(basename(f), "_")) for f in ARGS]
    allkeys = sort(unique(vcat(collect.(keys.(R))...)))
    for (title, col, fmt) in (("time [s]", k -> k == "inverse" ? "t_inv1" : "t_solve", sec),
                              ("peak memory [GB]", k -> k == "inverse" ? "dpeak_inv1" : "dpeak_solve", gb))
        println("\n", title); @printf("%-11s %3s %6s", "case", "n", "N"); foreach(n -> @printf(" %15s", n), names); println()
        for key in allkeys
            N = first(r[key]["N"] for r in R if haskey(r, key))
            @printf("%-11s %3d %6d", key[1], key[2], N)
            foreach(r -> @printf(" %15s", haskey(r, key) ? fmt(float(r[key][col(key[1])])) : "-"), R); println()
        end
    end
    println("\nmax |dE| of the lowest three energies against ", names[end])
    for key in allkeys
        haskey(R[end], key) || continue
        @printf("%-11s %3d       ", key[1], key[2])
        for r in R
            dE = haskey(r, key) ? maximum(abs(complex(r[key]["E$(i)_re"], r[key]["E$(i)_im"]) -
                 complex(R[end][key]["E$(i)_re"], R[end][key]["E$(i)_im"])) for i in 1:3) : NaN
            @printf(" %15s", isfinite(dE) ? @sprintf("%.1e", dE) : "-")
        end
        println()
    end
    exit()
end

A, B = load(ARGS[1]), load(ARGS[2])
println("old: ", basename(ARGS[1]), "\nnew: ", basename(ARGS[2]), "\n(ratio = new/old)\n")
@printf("%-11s %3s %6s | %8s %8s %6s | %8s %8s %6s | %7s %7s %6s | %8s\n", "case", "n", "N",
        "t_solve", "", "ratio", "t_build", "", "ratio", "dpeak", "[GB]", "ratio", "max|dE|")
for key in sort(collect(intersect(keys(A), keys(B))))
    a, b = A[key], B[key]
    inv = key[1] == "inverse"
    ts, tb, pk = inv ? ("t_inv1", "t_build", "dpeak_inv1") : ("t_solve", "t_build", "dpeak_solve")
    dE = maximum(abs(complex(b["E$(i)_re"], b["E$(i)_im"]) - complex(a["E$(i)_re"], a["E$(i)_im"])) for i in 1:3)
    @printf("%-11s %3d %6d | %s %s %s | %s %s %s | %s %s %s | %8.1e\n", key[1], key[2], b["N"],
            sec(a[ts]), sec(b[ts]), ratio(a[ts], b[ts]), sec(a[tb]), sec(b[tb]), ratio(a[tb], b[tb]),
            gb(a[pk]), gb(b[pk]), ratio(a[pk], b[pk]), dE)
    if inv
        @printf("%-11s %3d %6d | 6 energies: %s -> %s s, with inverse_prepare: %s s\n", "", key[2], b["N"],
                sec(a["t_inv6"]), sec(b["t_inv6"]), sec(b["t_inv6_prep"]))
    end
end

println("\npeak memory per N^2 (fit c0 + c1*N^2 over the sizes of a case) [B/N^2]")
for name in sort(unique(first.(collect(intersect(keys(A), keys(B))))))
    name == "inverse" && continue
    c1 = map((A, B)) do R
        ks = sort([k for k in keys(R) if k[1] == name])
        N = [R[k]["N"] for k in ks]
        (fit_c1(N, [float(R[k]["dpeak_solve"]) for k in ks]), fit_c1(N, [float(R[k]["dpeak_build"]) for k in ks]))
    end
    @printf("%-11s solve: %6.1f -> %6.1f   matrix build: %6.1f -> %6.1f\n", name, c1[1][1], c1[2][1], c1[1][2], c1[2][2])
end
