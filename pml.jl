function create_sfactor(wrange, s, omega, Nw, Nw_pml; m=3.5, lnR=-12)
    dw = (wrange[2]-wrange[1])/Nw; 
    Tw = Nw_pml*dw; 

    sig_max = -(m+1)*lnR/(2*eta0*Tw); 
    sig_w = l -> sig_max*(l/Tw)^m; 
    S = l -> 1-1im*sig_w(l)/(omega*epsilon0); 

    sfactor_array = ones(Complex128, Nw);

    for i in 1:Nw
        if s == "f"
            if i <= Nw_pml
                sfactor_array[i] = S(dw*(Nw_pml-i+0.5));
            elseif i > Nw - Nw_pml
                sfactor_array[i] = S(dw*(i-(Nw-Nw_pml)-0.5));
            end
        elseif s == "b"
            if i <= Nw_pml
                sfactor_array[i] = S(dw*(Nw_pml-i+1));
            elseif i > Nw - Nw_pml
                sfactor_array[i] = S(dw*(i-(Nw-Nw_pml)-1));
            end
        end
    end

    return sfactor_array
end

function S_create(omega, N, Npml, xrange, yrange)
    # Create the sfactor in each direction and for 'f' and 'b'
    s_vector_x_f = create_sfactor(xrange, "f", omega, N[1], Npml[1])
    s_vector_x_b = create_sfactor(xrange, "b", omega, N[1], Npml[1])
    s_vector_y_f = create_sfactor(yrange, "f", omega, N[2], Npml[2])
    s_vector_y_b = create_sfactor(yrange, "b", omega, N[2], Npml[2])

    # Fill the 2D space with layers of appropriate s-factors
    Sx_f_2D = zeros(Complex128, N);
    Sx_b_2D = zeros(Complex128, N);
    Sy_f_2D = zeros(Complex128, N);
    Sy_b_2D = zeros(Complex128, N);

    for i = 1:N[1]
        Sy_f_2D[i, :] = s_vector_y_f.^-1; 
        Sy_b_2D[i, :] = s_vector_y_b.^-1; 
    end

    for i = 1:N[2]
        Sx_f_2D[:, i] = s_vector_x_f.^-1;  
        Sx_b_2D[:, i] = s_vector_x_b.^-1; 
    end

    # Construct the 1D total s-array into a diagonal matrix
    Sx_f = spdiagm(Sx_f_2D[:]);
    Sx_b = spdiagm(Sx_b_2D[:]);
    Sy_f = spdiagm(Sy_f_2D[:]);
    Sy_b = spdiagm(Sy_b_2D[:]);

    return (Sx_f, Sx_b, Sy_f, Sy_b)
end