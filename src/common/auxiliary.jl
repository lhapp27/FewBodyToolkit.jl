# function often used in the examples to print a comparison of two arrays with differences and labels s1,s2
function comparison(num_arr,ref_arr,simax;s1="Numerical", s2="Reference", indexlist=1:simax)
    @printf("%-7s %-15s %-15s %-15s\n", "Index",  s1, s2, "Difference")
    for i in indexlist
        @printf("%-7d %-15.6f %-15.6f %-15.6f\n", i, num_arr[i], ref_arr[i], ref_arr[i] - num_arr[i])
    end
end;