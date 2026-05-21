# function to check the inputs for the three-body programs

function sanity_checks3B(phys_params)
    (;masses,species,interactions,parity) = phys_params
    
    error_code = 1 # default error code, if no errors are found, it will be set to 0
    if (lastindex(masses) !=3) || lastindex(species) !=3
        error("masses and/or species have wrong size, must be 3")
    else
        error_code = 0
    end
    
    fasb = findall(species .== :b)
    fasf = findall(species .== :f)
    if (lastindex(fasb) in [0,2,3]) && (lastindex(fasf) in [0,2,3])
        error_code = 0
    elseif lastindex(fasb) == 1 || lastindex(fasf) == 1
        error("Problem with symmetrization. Impossible number of identical particles. Must be either 0,2, or 3")
    end
    
    # tests required: ok, even if fasb or fasf is empty?
    if allequal(masses[fasb]) == false || allequal(masses[fasf]) == false
        error("Problem with symmetrization: species does not fit to m_arr")
    end
    
    return error_code
end
