# function to check the inputs for the three-body programs

function sanity_checks3B(phys_params)
    (;mass_arr,svals,vint_arr,parity) = phys_params
    
    if (lastindex(mass_arr) !=3) || lastindex(svals) !=3
        println("mass_arr and/or svals have wrong size, must be 3")
        error_code = 1
        return error_code
    else
        error_code = 0
    end
    
    fasb = findall(svals .== "b")
    fasf = findall(svals .== "f")
    if (lastindex(fasb) in [0,2,3]) && (lastindex(fasf) in [0,2,3])
        error_code = 0
    elseif lastindex(fasb) == 1 || lastindex(fasf) == 1
        println("Problem with symmetrization. Impossible number of identical particles. Must be either 0,2, or 3")
        error_code = 2
        return error_code
    end
    
    # tests required: ok, even if fasb or fasf is empty?
    if allequal(mass_arr[fasb]) == false || allequal(mass_arr[fasf]) == false
        println("Problem with symmetrization: svals does not fit to m_arr")
        error_code=3
        return error_code
    end
    
    return error_code
end
