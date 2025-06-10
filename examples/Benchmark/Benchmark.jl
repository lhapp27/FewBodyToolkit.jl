# Benchmark test for FewBodyToolkit vs FewBodyPhysics

using FewBodyToolkit, FewBodyPhysics, LinearAlgebra, BenchmarkTools

# 1. Example using FewBodyPhysics.jl:

function FBPrun(;output=true)
    masses = [1.0, 1.0, 1.0]
    psys = ParticleSystem(masses)
    
    K = Diagonal([1/2, 1/2, 1/2])
    K_transformed = psys.J * K * psys.J'
    
    w_list = [[1, -1, 0], [1, 0, -1], [0, 1, -1]]
    w_raw = [psys.U' * w for w in w_list]
    
    coeffs = [+1.0, -1.0, -1.0]
    
    let
        n_basis = 100
        b1 = default_b0(psys.scale)
        method = :psudorandom
        basis_fns = GaussianBase[]
        E₀_list = Float64[]
        
        for i in 1:n_basis
            bij = generate_bij(method, i, length(w_raw), b1)
            A = generate_A_matrix(bij, w_raw)
            push!(basis_fns, Rank0Gaussian(A))
            
            basis = BasisSet(basis_fns)
            ops = Operator[
            KineticEnergy(K_transformed);
            (CoulombPotential(c, w) for (c, w) in zip(coeffs, w_raw))...
            ]
            
            H = build_hamiltonian_matrix(basis, ops)
            S = build_overlap_matrix(basis)
            vals, vecs = solve_generalized_eigenproblem(H, S)
            #global c₀ = vecs[:, 1]
            E₀ = minimum(vals)
            
            push!(E₀_list, E₀)
            #println("Step $i: E₀ = $E₀")
        end
        
        if output==true
            E₀ = minimum(E₀_list)
            Eᵗʰ = -0.2620050702328
            ΔE = abs(E₀ - Eᵗʰ)
            println("FewBodyPhysics: ΔE = $ΔE")
        end
        
    end
    
end



# 2. Example using FewBodyToolkit.jl:

function FBTKrun(;output=true)
    # physical parameters:
    vee(r) = 1/r
    vep(r) = -1/r
    mass_arr = [1.0,1.0,1.0] # masses of particles: electron, electron, positron
    vint_arr = [[vep],[vep],[vee]] # v23,v31,v12
    phys_params = make_phys_params3B3D(;svals=["b","b","z"],vint_arr,spin_arr=[0,0,0],parity=+1,J_tot=0)
    
    # numerical parameters:
    gp = (;nmax=10,Nmax=10,r1=0.2,rnmax=20.0,R1=0.2,RNmax=20.0) # Tuple for GEM-Parameters
    lmax = 0;Lmax = 0; # maximum partial wave for r- and R- variables
    lmin = 0;Lmin = 0; # minimum partial wave for r- and R- variables
    num_params = make_num_params3B3D(;gem_params=gp,lmax=1,Lmax=1,kmax_interpol=100)
    
    energies_arr = ISGL_solve(phys_params,num_params);
    
    if output==true
        @show(energies_arr[1:4])
        E_ground = energies_arr[1]
        E_exact = -0.2620050702328
        ΔE = abs(E_ground - E_exact)
        println("FewBodyToolkit: ΔE = $ΔE")
    end
end

FBPrun()
FBTKrun()

println("\nBenchmarking FewBodyPhysics.jl vs FewBodyToolkit.jl:")
print("FewBodyPhysics.jl: ")
@time FBPrun(output=false)
print("FewBodyToolkit.jl: ")
@time FBTKrun(output=false)

println("")
#@benchmark FBTKrun(output=false)