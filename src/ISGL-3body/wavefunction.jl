## Function to calculate the three-body wave function in position space from the eigenvector obtained by the ISGL program

function wavefun(r,R,js,wavefun_arr,bf_info_arr,abvals_arr,box_size_arr,jmat)

    psi = 0.0

    index = 0;
    for (ii,val) in enumerate(abvals_arr) # for all boxes
        for jj = 1:box_size_arr[ii] # for all independent basis functions within a box
            index += 1
            for aa in val # for all ab-values (Jacobi-sets) within a box
                jm = jmat[aa,js]
                (la,La,nua,NUa,norm2) = bf_info_arr[index]

                contri = wavefun_arr[index]*norm2*(jm[1,1]*r + jm[1,2]*R)^la*(jm[2,1]*r + jm[2,2]*R)^La*exp(-nua*(jm[1,1]*r + jm[1,2]*R)^2 - NUa*(jm[2,1]*r + jm[2,2]*R)^2)
                @show([wavefun_arr[index],norm2,jm,la,La,nua,NUa])
                @show(contri)

                psi += contri
            end
        end
    end

    if index != lastindex(wavefun_arr)
        @show([index,lastindex(wavefun_arr)])
        error("Error in calculation of wave function. Index error.")
    end

    return psi
end
