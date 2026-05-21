# Direct low-level tests for GEM2B matrix-element dispatch and branch coverage

@testset "GEM2B Matrix Elements" begin
    # Small deterministic setup for direct matrix-element calls
    gamma_dict = Dict(
        0.5 => sqrt(pi),
        1.0 => 1.0,
        1.5 => 0.5 * sqrt(pi),
        2.0 => 1.0,
        2.5 => 0.75 * sqrt(pi),
        3.0 => 2.0,
    )

    buf_real = FewBodyToolkit.QuadGK.alloc_segbuf(Float64, Float64, Float64)
    buf_cplx = FewBodyToolkit.QuadGK.alloc_segbuf(Float64, ComplexF64, Float64)

    @testset "ContactPotential1D element_V branches" begin
        cp0 = ContactPotential1D(-1.0, 0.0)

        v_l0 = GEM2B.element_V(cp0, 0, 0.7, 1.1, gamma_dict, buf_real, 1, (-Inf, 0.0, Inf), 0.5)
        v_l1 = GEM2B.element_V(cp0, 1, 0.7, 1.1, gamma_dict, buf_real, 1, (-Inf, 0.0, Inf), 0.5)
        @test v_l0 < 0
        @test isapprox(v_l1, 0.0; atol=1e-14)

        cpz = ContactPotential1D(-1.0, 0.2)
        v_z = GEM2B.element_V(cpz, 0, 0.7, 1.1, gamma_dict, buf_real, 1, (-Inf, 0.0, Inf), 0.5)
        @test v_z < 0
    end

    @testset "CentralPotential and Function dispatch equivalence" begin
        vfun(r) = -1 / (1 + r^2)
        cp = CentralPotential(vfun)

        v1 = GEM2B.element_V(cp, 0, 0.5, 1.0, gamma_dict, buf_real, 3, (0.0, Inf), 1.0)
        v2 = GEM2B.element_V(vfun, 0, 0.5, 1.0, gamma_dict, buf_real, 3, (0.0, Inf), 1.0)
        @test isapprox(v1, v2; rtol=1e-10, atol=1e-12)
    end

    @testset "Invalid dimension branches" begin
        @test_throws ErrorException GEM2B.element_T(0, 0.7, 1.1, 1.0, 1.0, 4)

        nu_arr = [0.4, 0.9]
        V = zeros(2, 2)
        @test_throws ErrorException GEM2B.MatrixV(V, 0, nu_arr, [GaussianPotential(-1.0, 0.5)], gamma_dict, buf_real, false, 0.0, false, 4)

        WD = zeros(2, 2)
        @test_throws ErrorException GEM2B.MatrixWD(WD, 0, nu_arr, r -> 0.0, [1, r -> 0.0], gamma_dict, buf_real, false, 0.0, false, 4)
    end

    @testset "MatrixT and MatrixV full-fill path" begin
        nmax = 3
        nu_cr = zeros(ComplexF64, 2 * nmax)
        GEM2B.buildnu(nu_cr, 0.3, 3.0, nmax)
        nu_cr .*= (1 + 0.4im)
        nu_cr[nmax + 1:2 * nmax] .= conj.(nu_cr[1:nmax])

        T = zeros(ComplexF64, 2 * nmax, 2 * nmax)
        GEM2B.MatrixT(T, 0, nu_cr, 1.0, 1.0, true, 5.0, true, 3)
        @test any(abs(T[i, j]) > 0 for i in 1:size(T, 1) for j in i + 1:size(T, 2))

        V = zeros(ComplexF64, 2 * nmax, 2 * nmax)
        GEM2B.MatrixV(V, 0, nu_cr, [GaussianPotential(-1.0, 0.5)], gamma_dict, buf_cplx, true, 5.0, true, 3)
        @test any(abs(V[i, j]) > 0 for i in 1:size(V, 1) for j in i + 1:size(V, 2))
    end

    @testset "MatrixWD derivative branches" begin
        nmax = 3
        nu_arr = zeros(nmax)
        GEM2B.buildnu(nu_arr, 0.4, 2.0, nmax)

        WCC(r) = 0.05 * exp(-r^2)
        pref(r) = 0.1 * exp(-r^2)

        for n in 1:4
            WDn = zeros(nmax, nmax)
            GEM2B.MatrixWD(WDn, 0, nu_arr, WCC, [n, pref], gamma_dict, buf_real, false, 0.0, true, 3)
            @test all(isfinite.(WDn))
        end

        WD_bad = zeros(nmax, nmax)
        @test_throws ErrorException GEM2B.MatrixWD(WD_bad, 0, nu_arr, WCC, [5, pref], gamma_dict, buf_real, false, 0.0, true, 3)
    end
end
