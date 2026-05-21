# function to check the inputs for the ISGL program of sanity:

function sanity_checks(phys_params)
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
    
    # tests required: ok, even if fasb or fasf is empty?
    if allequal(masses[fasb]) == false || allequal(masses[fasf]) == false
        println("Problem with symmetrization: species does not fit to m_arr")
        error_code=3
        return error_code
    end
    
    return error_code
end
