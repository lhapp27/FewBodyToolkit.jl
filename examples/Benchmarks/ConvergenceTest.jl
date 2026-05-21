## Convergence check for the 2+1 example code.

using Printf, FewBodyToolkit

# Same helper as the example
function exfun(mr)
    if mr == 2.2
        return [-2.1966, -1.0520]
    elseif mr == 12.4
        return [-2.5963, -1.4818, -1.1970, -1.0377, -1.0002]
    elseif mr == 22.2
        return [-2.7515, -1.6904, -1.3604, -1.1479, -1.0525, -1.0040]
    else
        error("Unknown mass ratio: $mr")
    end
end

# -------- Fixed setup (copied from the example logic) --------
massratio = 22.2
masses = [1.0, massratio, massratio]
mur = 1 / (1 / masses[1] + 1 / masses[2])

v0 = -1.0
mu_g = 1.0
vg = GaussianPotential(v0, mu_g)

# 2B inverse problem for Gaussian interaction
phys_params2B = make_phys_params2B(; mur, interactions=[vg], dim=1)
num_params2B = make_num_params2B(; gem_params=(; nmax=16, r1=1.0, rnmax=120.0))
stateindex = 1
target_e2 = -1e-3

pps, nps, vscale = GEM2B.v0GEMOptim(phys_params2B, num_params2B, stateindex, target_e2)
vgscaled = GaussianPotential(v0 * vscale, mu_g)
pps = make_phys_params2B(; mur, interactions=[vgscaled], dim=1)
e2s = GEM2B.GEM2B_solve(pps, nps)
e2g = abs(e2s[1])

# 2B contact interaction (as in the example)
vc = ContactPotential1D(-sqrt(-2 * target_e2), 0.0)
ppc = make_phys_params2B(; mur, interactions=[vc], dim=1)
npc = make_num_params2B(; gem_params=(; nmax=16, r1=1.0, rnmax=120.0))
r1cs, rnmaxcs, e2copt = GEM_Optim_2B(ppc, npc, stateindex)
e2c = abs(e2copt)

function solve_all_eps(n; masses, vgscaled, vc, e2g, e2c, r1g, rnmaxg, r1c, rnmaxc)
    # Bosons, Gaussian
    interactions = [[], [vgscaled], [vgscaled]]
    pp3B = make_phys_params3B1D(; masses=masses, species=[:x,:b,:b], interactions=interactions)
    np3B = make_num_params3B1D(; gem_params=(; nmax=n, r1=r1g, rnmax=rnmaxg, Nmax=n, R1=1.5, RNmax=250.0))
    e3 = GEM3B1D.GEM3B1D_solve(pp3B, np3B)
    epsilon = e3 / e2g

    # Bosons, Contact
    vint_arrC = [[], [vc], [vc]]
    pp3BC = make_phys_params3B1D(; masses=masses, species=[:x,:b,:b], interactions=vint_arrC)
    np3BC = make_num_params3B1D(; gem_params=(; nmax=n, r1=r1c, rnmax=rnmaxc, Nmax=n, R1=1.5, RNmax=250.0))
    e3c = GEM3B1D.GEM3B1D_solve(pp3BC, np3BC)
    epsilonC = e3c / e2c

    # Fermions, Gaussian
    pp3BF = make_phys_params3B1D(; masses=masses, species=[:x,:f,:f], interactions=interactions, parity=-1)
    np3BF = make_num_params3B1D(;
        gem_params=(; nmax=n, r1=r1g, rnmax=rnmaxg, Nmax=n, R1=1.5, RNmax=250.0),
        lmin=0, Lmin=0, lmax=1, Lmax=1
    )
    e3_F = GEM3B1D.GEM3B1D_solve(pp3BF, np3BF)
    epsilon_F = e3_F / e2g

    # Fermions, Contact
    pp3BCF = make_phys_params3B1D(; masses=masses, species=[:x,:f,:f], interactions=vint_arrC, parity=-1)
    np3BCF = make_num_params3B1D(;
        gem_params=(; nmax=n, r1=r1c, rnmax=rnmaxc, Nmax=n, R1=1.5, RNmax=250.0),
        lmin=0, Lmin=0, lmax=1, Lmax=1
    )
    e3c_F = GEM3B1D.GEM3B1D_solve(pp3BCF, np3BCF)
    epsilonC_F = e3c_F / e2c

    return (epsilon=epsilon, epsilonC=epsilonC, epsilon_F=epsilon_F, epsilonC_F=epsilonC_F)
end

first_or_nan(v) = isempty(v) ? NaN : v[1]

# Max relative change on common bound-state subset between consecutive n
function max_rel_change(curr::Vector{<:Real}, prev::Vector{<:Real})
    m = min(lastindex(curr), lastindex(prev), 3)
    m == 0 && return NaN
    den = max.(abs.(curr[1:m]), 1e-12)
    return maximum(abs.((curr[1:m] .- prev[1:m]) ./ den))
end

nvals = [6, 8, 10, 12, 14, 16]

println("Convergence test with nmax = Nmax")
println("massratio = ", massratio)
println("")
println("Columns:")
println(" n |  eps[1]      epsC[1]     epsF[1]     epsCF[1]   |  d_eps      d_epsC     d_epsF     d_epsCF")
println("-"^100)

let prev = nothing
    for n in nvals
        r = solve_all_eps(
            n;
            masses=masses, vgscaled=vgscaled, vc=vc, e2g=e2g, e2c=e2c,
            r1g=nps.gem_params.r1, rnmaxg=nps.gem_params.rnmax,
            r1c=r1cs, rnmaxc=rnmaxcs
        )

        if prev === nothing
            @printf("%2d | % .8f % .8f % .8f % .8f | %8s %8s %8s %8s\n",
                n,
                first_or_nan(r.epsilon), first_or_nan(r.epsilonC),
                first_or_nan(r.epsilon_F), first_or_nan(r.epsilonC_F),
                "-", "-", "-", "-"
            )
        else
            d1 = max_rel_change(r.epsilon, prev.epsilon)
            d2 = max_rel_change(r.epsilonC, prev.epsilonC)
            d3 = max_rel_change(r.epsilon_F, prev.epsilon_F)
            d4 = max_rel_change(r.epsilonC_F, prev.epsilonC_F)

            @printf("%2d | % .8f % .8f % .8f % .8f | % .3e % .3e % .3e % .3e\n",
                n,
                first_or_nan(r.epsilon), first_or_nan(r.epsilonC),
                first_or_nan(r.epsilon_F), first_or_nan(r.epsilonC_F),
                d1, d2, d3, d4
            )
        end

        prev = r
    end
end

println("")
println("Reference ratios (contact limit) for this mass ratio:")
refB = exfun(massratio)[1:2:end]
refF = exfun(massratio)[2:2:end]
println("bosons:   ", refB)
println("fermions: ", refF)