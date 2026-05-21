# Tests for wavefunction utilities and preallocation type combinations

@testset "GEM2B Wavefunctions and Prealloc" begin
    # Reuse a small Coulomb setup for deterministic and fast checks
    v_coulomb(r) = -1.0 / r
    phys3 = make_phys_params2B(; interactions=[v_coulomb], dim=3)
    num3 = make_num_params2B(; gem_params=(nmax=8, r1=0.12, rnmax=30.0))

    energies, wfs = GEM2B.GEM2B_solve(phys3, num3; return_wavefunctions=true)
    @test all(isfinite.(real.(energies[1:4])))

    @testset "wavefun_point and wavefun_arr" begin
        rgrid = collect(range(0.1, 2.0; length=12))

        psi_grid = GEM2B.wavefun_arr(rgrid, phys3, num3, wfs[:, 1])
        @test length(psi_grid) == length(rgrid)
        @test all(isfinite.(real.(psi_grid)))

        # Build nu directly to call wavefun_point
        nu_arr = zeros(num3.gem_params.nmax)
        GEM2B.buildnu(nu_arr, num3.gem_params.r1, num3.gem_params.rnmax, num3.gem_params.nmax)
        psi_pt = GEM2B.wavefun_point(1.0, nu_arr, wfs[:, 1], phys3.lmax, phys3.dim)
        @test isfinite(real(psi_pt))
    end

    @testset "complex-ranged wavefun_arr" begin
        num_cr = make_num_params2B(; gem_params=(nmax=6, r1=0.2, rnmax=8.0), complex_range_freq=0.6)
        energies_cr, wfs_cr = GEM2B.GEM2B_solve(phys3, num_cr; return_wavefunctions=true, complex_ranged=true)
        @test all(isfinite.(real.(energies_cr[1:4])))

        rgrid = collect(range(0.1, 1.5; length=8))
        psi_grid_cr = GEM2B.wavefun_arr(rgrid, phys3, num_cr, wfs_cr[:, 1]; complex_ranged=true)
        @test length(psi_grid_cr) == length(rgrid)
        @test all(isfinite.(real.(psi_grid_cr)))
    end

    @testset "PreallocStruct2B combinations" begin
        np = make_num_params2B(; gem_params=(nmax=5, r1=0.2, rnmax=10.0))

        pa_rr = GEM2B.PreallocStruct2B(np, false, false)
        @test eltype(pa_rr.nu_arr) == Float64
        @test eltype(pa_rr.S) == Float64
        @test eltype(pa_rr.energies) == Float64
        @test size(pa_rr.wavefunctions) == (5, 5)

        pa_cr = GEM2B.PreallocStruct2B(np, true, false)
        @test eltype(pa_cr.nu_arr) == ComplexF64
        @test eltype(pa_cr.S) == ComplexF64
        @test eltype(pa_cr.energies) == Float64
        @test size(pa_cr.wavefunctions) == (10, 10)

        pa_csm = GEM2B.PreallocStruct2B(np, false, true)
        @test eltype(pa_csm.nu_arr) == ComplexF64
        @test eltype(pa_csm.S) == Float64
        @test eltype(pa_csm.energies) == ComplexF64
        @test size(pa_csm.wavefunctions) == (5, 5)

        pa_both = GEM2B.PreallocStruct2B(np, true, true)
        @test eltype(pa_both.nu_arr) == ComplexF64
        @test eltype(pa_both.S) == ComplexF64
        @test eltype(pa_both.energies) == ComplexF64
        @test size(pa_both.wavefunctions) == (10, 10)
    end
end
