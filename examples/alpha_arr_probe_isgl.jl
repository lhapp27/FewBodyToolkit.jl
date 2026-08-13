using FewBodyToolkit
using LinearAlgebra
using Printf
using StaticArrays

const ISGL = FewBodyToolkit.ISGL

"""
Build ISGL Jacobi transformation matrices exactly via the internal ISGL precompute helper.
"""
function build_jmat(masses::Vector{Float64})
    jmat = Matrix{SMatrix{2,2,Float64,4}}(undef, 3, 3)
    ISGL.precompute_jmat(jmat, masses)
    return jmat
end

@views function buildnu!(r1,rnmax,nmax,nu_arr,complex_ranged::Bool=false,complex_range_freq=0.0)
    nu_arr[1] = 1 /r1^2;
    nmax >1 && @. nu_arr[2:nmax] = 1/r1^2 * (r1/rnmax)^(2*((2:nmax)-1)/(nmax-1))
    if complex_ranged
        nu_arr[1:nmax] .*= (1 + complex_range_freq*im)
        nu_arr[nmax+1:2*nmax] .= conj.(nu_arr[1:nmax])
    end
    #return nu_arr # we dont need to return since we overwrite the array in-place
end

"""
Build Gaussian ranges with the same buildnu! routine used in ISGL precomputation.
"""
function build_ranges(r1, rnmax, nmax; complex_ranged::Bool=false, complex_range_freq::Float64=0.0)
    if complex_ranged
        nu = Vector{ComplexF64}(undef, 2 * nmax)
    else
        nu = Vector{Float64}(undef, nmax)
    end
    buildnu!(r1, rnmax, nmax, nu, complex_ranged, complex_range_freq)
    return nu
end

"""
Convert alpha bounds (tempmin/tempmax in the interpolation routine sense) to the alpha grid.
"""
function alpha_grid_from_bounds(tempmin::Float64, tempmax::Float64, kmax::Int)
    tempmin > 0 || error("tempmin must be positive, got $tempmin")
    tempmax > 0 || error("tempmax must be positive, got $tempmax")
    tempmax >= tempmin || error("Expected tempmax >= tempmin, got tempmin=$tempmin tempmax=$tempmax")
    
    min_r_eff = 1 / sqrt(tempmax)
    max_r_eff = 1 / sqrt(tempmin)
    alpha_arr = zeros(kmax)
    buildnu!(max_r_eff, min_r_eff, kmax, alpha_arr)
    return alpha_arr
end

"""
Current ISGL full scan logic (parameterized by the scan projection used for nu/NU).
projection = :real  -> use real.(nu), real.(NU) (current implementation)
projection = :abs   -> use abs.(nu),  abs.(NU)
"""
function alpha_bounds_full(jmat, nu_arr, NU_arr; projection::Symbol=:real)
    nu_scan = projection === :real ? real.(nu_arr) : abs.(nu_arr)
    NU_scan = projection === :real ? real.(NU_arr) : abs.(NU_arr)
    
    first_iter = true
    tempmin = 0.0
    tempmax = 0.0
    
    for a in 1:3
        for b in 1:3
            for c in 1:3
                for nua in nu_scan
                    for NUa in NU_scan
                        for nub in nu_scan
                            for NUb in NU_scan
                                Aa = SA[nua 0.0; 0.0 NUa]
                                Ab = SA[nub 0.0; 0.0 NUb]
                                tempA = transpose(jmat[a, c]) * Aa * jmat[a, c] + transpose(jmat[b, c]) * Ab * jmat[b, c]
                                temp = real(det(tempA) / tempA[2, 2])
                                
                                if first_iter
                                    tempmin = temp
                                    tempmax = temp
                                    first_iter = false
                                else
                                    temp < tempmin && (tempmin = temp)
                                    temp > tempmax && (tempmax = temp)
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    # same numerical padding style as ISGL code
    tempmin -= 100 * eps()
    tempmax += 100 * eps()
    
    return tempmin, tempmax
end

"""
1D-style factorized bounds adapted to ISGL:
- min/max over pure j-matrix geometry
- scale by global min/max Gaussian widths derived from endpoints

For complex ranges we offer two endpoint projections:
endpoint_projection = :real or :abs.
"""
function alpha_bounds_factorized(jmat, r1, rnmax, R1, RNmax; endpoint_projection::Symbol=:real)
    first_iter = true
    tempmin2 = 0.0
    tempmax2 = 0.0
    
    for a in 1:3
        for b in 1:3
            for c in 1:3
                tempA2 = transpose(jmat[a, c]) * jmat[a, c] + transpose(jmat[b, c]) * jmat[b, c]
                temp2 = det(tempA2) / tempA2[2, 2]
                
                if first_iter
                    tempmin2 = temp2
                    tempmax2 = temp2
                    first_iter = false
                else
                    temp2 < tempmin2 && (tempmin2 = temp2)
                    temp2 > tempmax2 && (tempmax2 = temp2)
                end
            end
        end
    end
    
    rmin = endpoint_projection === :real ? min(real(r1), real(R1)) : min(abs(r1), abs(R1))
    rmax = endpoint_projection === :real ? max(real(rnmax), real(RNmax)) : max(abs(rnmax), abs(RNmax))
    
    tempmin = tempmin2 / rmax^2
    tempmax = tempmax2 / rmin^2
    
    # same safety padding
    tempmin -= 100 * eps()
    tempmax += 100 * eps()
    
    return tempmin, tempmax
end

function alpha_full(jmat, nu_arr, NU_arr)
    
    alpha_arr=zeros(ComplexF64,1)
    alpha_arr2=zeros(ComplexF64,1)
    
    for a in 1:3
        for b in 1:3
            for c in 1:3
                for nua in nu_arr
                    for NUa in NU_arr
                        for nub in nu_arr
                            for NUb in NU_arr
                                Aa = SA[nua 0.0; 0.0 NUa]
                                Ab = SA[nub 0.0; 0.0 NUb]
                                
                                tempA = transpose(jmat[a, c]) * Aa * jmat[a, c] + transpose(jmat[b, c]) * Ab * jmat[b, c]
                                temp = det(tempA) / tempA[2, 2]
                                push!(alpha_arr, temp)
                                
                                tempA2 = transpose(jmat[a, c]) * conj.(Aa) * jmat[a, c] + transpose(jmat[b, c]) * Ab * jmat[b, c]
                                temp2 = det(tempA2) / tempA2[2, 2]
                                push!(alpha_arr2, temp2)
                                
                            end
                        end
                    end
                end
            end
        end
    end
    
    return alpha_arr, alpha_arr2
end

function summarize_case(title, jmat, r1, rnmax, R1, RNmax, nu_arr, NU_arr; kmax=16, full_projection::Symbol=:real, factor_projection::Symbol=:real)
    @printf("\n=== %s ===\n", title)
    
    tmin_full, tmax_full = alpha_bounds_full(jmat, nu_arr, NU_arr; projection=full_projection)
    tmin_fac, tmax_fac = alpha_bounds_factorized(jmat, r1, rnmax, R1, RNmax; endpoint_projection=factor_projection)
    
    alpha_full = alpha_grid_from_bounds(tmin_full, tmax_full, kmax)
    alpha_fac = alpha_grid_from_bounds(tmin_fac, tmax_fac, kmax)
    
    @printf("full bounds      : [%.8e, %.8e]\n", tmin_full, tmax_full)
    @printf("factorized bounds: [%.8e, %.8e]\n", tmin_fac, tmax_fac)
    
    conservative = (tmin_fac <= tmin_full) && (tmax_fac >= tmax_full)
    @printf("factorized brackets full? %s\n", conservative ? "YES" : "NO")
    
    @printf("alpha_full first/last: %.8e  %.8e\n", first(alpha_full), last(alpha_full))
    @printf("alpha_fac  first/last: %.8e  %.8e\n", first(alpha_fac), last(alpha_fac))
    
    @printf("relative diff min bound: %.4e\n", abs(tmin_fac - tmin_full) / abs(tmin_full))
    @printf("relative diff max bound: %.4e\n", abs(tmax_fac - tmax_full) / abs(tmax_full))
end

function run_probe(; nmax=4, Nmax=4, kmax=20)
    masses = [1.0, 2000.0, 2000.0]
    jmat = build_jmat(masses)
    
    # Base real endpoints
    r1 = 0.1
    rnmax = 5.0
    R1 = 0.5
    RNmax = 10.0
    
    # 1) Purely real ranges
    nu_real = build_ranges(r1, rnmax, nmax; complex_ranged=false)
    NU_real = build_ranges(R1, RNmax, Nmax; complex_ranged=false)
    summarize_case(
    "REAL ranges: full(real-scan) vs factorized(real-endpoints)",
    jmat, r1, rnmax, R1, RNmax, nu_real, NU_real;
    kmax=kmax, full_projection=:real, factor_projection=:real,
    )
    
    # 2) Complex-ranged inputs, compared with current ISGL-style real projection
    freq = 0.9
    nu_cr = build_ranges(r1, rnmax, nmax; complex_ranged=true, complex_range_freq=freq)
    NU_cr = build_ranges(R1, RNmax, Nmax; complex_ranged=true, complex_range_freq=freq)
    r1_cr = r1 / sqrt(1 + freq^2)
    rnmax_cr = rnmax / sqrt(1 + freq^2)
    R1_cr = R1 / sqrt(1 + freq^2)
    RNmax_cr = RNmax / sqrt(1 + freq^2)
    
    summarize_case(
    "COMPLEX-ranged (freq=$freq): full(real-projection) vs factorized(real-endpoints)",
    jmat, r1, rnmax, R1, RNmax, nu_cr, NU_cr;
    kmax=kmax, full_projection=:real, factor_projection=:real,
    )
    
    summarize_case(
    "COMPLEX-ranged (freq=$freq): full(abs-projection) vs factorized(abs-endpoints)",
    jmat, r1_cr, rnmax_cr, R1_cr, RNmax_cr, nu_cr, NU_cr;
    kmax=kmax, full_projection=:abs, factor_projection=:abs,
    )
    
    println("\nDone. If desired, we can next plug in your exact production masses and gem_params.")
end

run_probe()

println("1. factorized (i.e. simpler approach) is conservative, but probably fine? works even for very large mass ratios and r1,rnmax, R1, RNmax discrepancies.")
println("2. why the hell should alpha be real?! i think it takes on complex values?!... ah man, i think we need to debug the 1D code a bit more to make sure about the actual values of alpha, and the interpolated ones! maybe interpolation actually works, even though it was just defined on the real axis?! but that would be a bit scary to rely on!")
