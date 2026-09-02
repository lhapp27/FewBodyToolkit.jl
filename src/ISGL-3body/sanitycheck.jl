# function to check the inputs for the ISGL program of sanity:

function sanity_checks(phys_params,complex_ranged::Bool=false,observ_params=nothing)
    (;masses,species,interactions,J_tot,parity) = phys_params

    if (lastindex(masses) !=3) || lastindex(species) !=3
        println("masses and/or species have wrong size, must be 3")
        error_code = 1
        return error_code
    else
        error_code = 0
    end
    
    fasb = findall(species .== :b)
    fasf = findall(species .== :f)
    if (lastindex(fasb) in [0,2,3]) && (lastindex(fasf) in [0,2,3])
        error_code = 0
    elseif lastindex(fasb) == 1 || lastindex(fasf) == 1
        println("Problem with symmetrization. Impossible number of identical particles. Must be either 0,2, or 3")
        error_code = 2
        return error_code
    end
    
    # power-law interactions: the radial integral int r^(2n+2) V(r) exp(-alpha r^2) dr
    # with n=0 converges only for p > -3.
    for c in 1:lastindex(interactions)
        for vint in interactions[c]
            if vint isa PowerLawPotential && vint.p <= -3.0
                println("PowerLawPotential: p = $(vint.p) is too singular for ISGL; requires p > -3")
                error_code = 4
                return error_code
            end
        end
    end
    
    # complex-ranged basis functions:
    # The generic path for a central potential obtains its radial integral numerically and then
    # interpolates at etaprc, which is complex once the ranges are. That is covered: the mesh is a
    # sector (see theta_mesh), so any non-singular central potential works. Observables, however,
    # still use the plain one-dimensional mesh in log(alpha) and remain unsupported.
    if complex_ranged
        if !isnothing(observ_params)
            (;centobs_arr,R2_arr) = observ_params
            if any(.!isempty.(centobs_arr)) || any(R2_arr .!= 0)
                println("Observables (centobs_arr, R2_arr) are not supported together with complex-ranged basis functions.")
                error_code = 6
                return error_code
            end
        end
    end

    # tests required: ok, even if fasb or fasf is empty?
    if allequal(masses[fasb]) == false || allequal(masses[fasf]) == false
        println("Problem with symmetrization: species does not fit to m_arr")
        error_code=3
        return error_code
    end
    
    return error_code
end
