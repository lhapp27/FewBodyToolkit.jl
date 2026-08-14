using CairoMakie

include("PerformanceScaling_nthreads=1.jl")

# extract columns
nmax_r, time_r, gc_r, mem_r, allocs_r, eig_r = eachcol(results)
nmax_a, time_a, gc_a, mem_a, allocs_a, eig_a = eachcol(results_analytical)

time_num = time_r ./ 1e9  # convert from ns to s
time_anal = time_a ./ 1e9  # convert from ns to s
memory_num = mem_r ./ 1e6  # convert from bytes to MB
memory_anal = mem_a ./ 1e6  # convert from bytes to MB

# Number of basis functions
nbasis = 3 .* nmax_r .^2

# Create figure
colors = Makie.wong_colors()
set_theme!(theme_latexfonts())
fig = Figure(fontsize=18,size=(800,400))

# Time plot
ax1 = Axis(fig[1,1],
    xlabel="Number of basis functions",
    ylabel="Time [s]",
    xscale=log2,yscale=log10,
    xminorticksvisible=true,yminorticksvisible=true,
    xminorticks = IntervalsBetween(5),yminorticks = IntervalsBetween(10),
    xgridvisible=false,ygridvisible=false,
    xtickformat = "{:d}",
    limits=(nothing,(10.0^-3,10.0^2)),
)
ax1top = Axis(fig[1, 1],
    xaxisposition = :top,
    xlabel = "n_max",
    xscale=log2,
    xminorticksvisible=true,
    xminorticks = IntervalsBetween(4),
    xgridvisible=false,ygridvisible=false,
    xtickformat = "{:d}",
    limits=(nothing,(1,2)),
)
hideydecorations!(ax1top)

scatterlines!(ax1, nbasis, time_num, color=:blue, markersize = 12, label="Numerical")
scatterlines!(ax1, nbasis, time_anal, color=:red, markersize = 12, label="Analytical",marker=:diamond)
scatter!(ax1top,nmax_a,fill(3,lastindex(nmax_a)))
axislegend(ax1, position=:rb)

# Memory plot
ax2 = Axis(fig[1,2],
    xlabel="Number of basis functions",
    ylabel= L"Memory [$10^6$ B]",
    xscale=log2,yscale=log10,
    xminorticksvisible=true,yminorticksvisible=true,
    xminorticks = IntervalsBetween(5),yminorticks = IntervalsBetween(10),
    xgridvisible=false,ygridvisible=false,
    xtickformat = "{:d}",
    limits=(nothing,(10.0^-1,10.0^4)),
)
ax2top = Axis(fig[1, 2],
    xaxisposition = :top,
    xlabel = "n_max",
    xscale=log2,
    xminorticksvisible=true,
    xminorticks = IntervalsBetween(4),
    xgridvisible=false,ygridvisible=false,
    xtickformat = "{:d}",
    limits=(nothing,(1,2)),
)
hideydecorations!(ax2top)

scatterlines!(ax2, nbasis, memory_num, color=:blue, markersize = 12, label="Numerical")
scatterlines!(ax2, nbasis, memory_anal, color=:red, markersize = 12, label="Analytical",marker=:diamond)
scatter!(ax2top,nmax_a,fill(3,lastindex(nmax_a)))
axislegend(ax2, position=:rb)

fig;
#save("bench.pdf",fig)