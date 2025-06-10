# Example script for using the complex-scaling method for resonance calculations in a two-body system. Comparison with https://doi.org/10.1007/s00601-023-01808-x

using FewBodyToolkit, Plots

# This function computes the complex energy spectrum for a given lambda of the potential defined int eh reference
function csm2b(lambda)
    function v(r)
        r > 50.0 && return 0.0
        return lambda*(678.1*exp(-2.55*r) - 166.0*exp(-0.68*r))/r
    end
    
    mur = 1/(2*27.647)
    pp = make_phys_params2B(;mur,vint_arr=[v],dim=3,lmin=1,lmax=1)
    np = make_num_params2B(;gem_params=(nmax=30, r1=0.3, rnmax=30.8),theta_csm=40.0)
    
    e2 = GEM2B_solve(pp,np,csm_bool=1)
    
    return e2
end

# This function filters the resonances of an array A in the complex energy plane, which are restricted by the triangle defined by
# - the location of the top left corner e21
# - the location position the top right corner e22
# - the angle between the horizontal cathetus e21-e22 and the hypothenuse from e21 to the bottom right corner e23 (location of e23 is computed within the function)
function is_in_triangle(A, e21, e22, angle)
    b_in = [false for _ in A]
    e23_real = real(e22)
    e23_imag = - abs(real(e22) - real(e21)) * tan(angle*pi/180)
    #println(e23_imag)
    
    for i = 1:lastindex(A)
        x = real(A[i])
        y = imag(A[i])
        
        if (x >= real(e21)) && (x <= real(e22)) && (y <= imag(e22)) && (y >= e23_imag) &&
            ((y - imag(e21)) >= -tan(angle*pi/180) * abs(x - real(e21)))
            b_in[i] = true
        end
        #println(x,",",y,",",e21,",",e22,",",e23_real + im*e23_imag)
    end
    
    return b_in
end

# Arrays and loop to store the results for various lambda values of the article.
lambda_arr = [1.0,1.25,1.5,1.75]
e2_arr = zeros(ComplexF64,30,lastindex(lambda_arr))
reso_arr = zeros(ComplexF64,lastindex(lambda_arr))
for (il,lambda) in enumerate(reverse(lambda_arr))
    e2_arr[:,il] = csm2b(lambda)

    reso = e2_arr[is_in_triangle(e2_arr[:,il], 0.0, 2.0, 2*40),il]
    if lambda == 1.75 # in this case there is only a bound state, no resonance
        reso_arr[il] = e2_arr[1,il]
    else
        reso_arr[il] = reso[1]
    end
end

reso_exact = [-1.7914+0.0*im,0.0933-0.0152*im,0.9712-0.7446*im,1.267-2.002*im,0.7429-3.530*im]

println("Real part:")
comparison(real.(reso_arr), real.(reso_exact),4)

println("\nImaginary part:")
comparison(imag.(reso_arr), imag.(reso_exact),4)
