# Script for testing the scaling of performance (time and memory) with increasing number of basis functions.
# Moreover, the performance between numerical and analytical treatment of the potential is compared

using FewBodyToolkit, BenchmarkTools, Dates, Printf

# function for the numerical approach
function perftest(pp,np)
    y = ISGL_solve(pp,np)[1];
    x = @benchmark ISGL_solve($pp,$np)
    t,g,m,a = mean(x).time,mean(x).gctime,mean(x).memory,mean(x).allocs
    nn = np.gem_params.nmax
    @show([nn,t,g,m,a,y])
    return (nn, t, g, m, a, y)
end

# function for the analytical approach
function perftest_analytical(pp,np)
    y = ISGL_solve(pp,np)[1];
    x = @benchmark ISGL_solve($pp,$np)
    t,g,m,a = mean(x).time,mean(x).gctime,mean(x).memory,mean(x).allocs
    nn = np.gem_params.nmax
    @show([nn,t,g,m,a,y])
    return (nn, t, g, m, a, y)
end

vg(r) = -10.0*exp(-r^2)
pp = make_phys_params3B3D(;mass_arr=[1.0,2.0,3.0],svals=["x","y","z"],vint_arr=[[vg],[vg],[vg]])

vga = GaussianPotential(-10.0,1.0)
ppa = make_phys_params3B3D(;mass_arr=[1.0,2.0,3.0],svals=["x","y","z"],vint_arr=[[vga],[vga],[vga]])

nns = 4:2:40
results = zeros(length(nns), 6)
results_analytical = zeros(length(nns), 6)
for (in,nn) in enumerate(nns)
    gp = (;nmax=nn,Nmax=nn,r1=0.1,rnmax=100.0,R1=0.1,RNmax=100.0)
    np = make_num_params3B3D(;gem_params=gp,kmax_interpol=2000,lmin=0,lmax=0,Lmin=0,Lmax=0,threshold=10^-10)

    results[in,:] .= perftest(pp,np)
    results_analytical[in,:] .= perftest_analytical(ppa,np)
end

nthreads = Threads.nthreads()

filename = "PerformanceScaling_nthreads=$(nthreads).jl"
open(filename, "w") do io
    println(io, "# date = $(Dates.now())")
    println(io, "# columns correspond to: nmax, time, gc_time, memory, allocs, eigenvalue")
    println(io, "results = [")
    for row in eachrow(results)
        print(io, "  ")
        for val in row
            @printf(io, "%11.5g  ", val)
        end
        println(io)
    end
    println(io, "]\n")

    println(io, "results_analytical = [")
    for row in eachrow(results_analytical)
        print(io, "  ")
        for val in row
            @printf(io, "%11.5g  ", val)
        end
        println(io)
    end
    println(io, "]")
end