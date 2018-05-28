function create_sfactor(w::Direction, s::DerivativeDirection, g::Grid, ω::Float ; m=3.5, lnR=-12)
    (ϵ₀, μ₀, c₀) = normalize_parameters(g);
    dw = dh(g, w);
    Tw = g.Npml[Int(w)]*dw;
    Nw = size(g, Int(w));
    Nw_pml = g.Npml[Int(w)];

    σmax = -(m+1)*lnR/(2*η₀*Tw);
    σw = l -> σmax*(l/Tw)^m;
    S = l -> 1-1im*σw(l)/(ω*ϵ₀);

    sfactor_array = ones(Complex128, Nw);

    for i in 1:Nw
        if s == Forward
            if i <= Nw_pml
                sfactor_array[i] = S(dw*(Nw_pml-i+0.5));
            elseif i > Nw - Nw_pml
                sfactor_array[i] = S(dw*(i-(Nw-Nw_pml)-0.5));
            end
        elseif s == Backward
            if i <= Nw_pml
                sfactor_array[i] = S(dw*(Nw_pml-i+1));
            elseif i > Nw - Nw_pml
                sfactor_array[i] = S(dw*(i-(Nw-Nw_pml)-1));
            end
        end
    end

    return sfactor_array
end

function S_create(g::Grid, ω::Float)
    # Create the sfactor in each direction and for forward and backward
    s_vector_x_f = create_sfactor(x̂, Forward,  g, ω);
    s_vector_x_b = create_sfactor(x̂, Backward, g, ω);
    s_vector_y_f = create_sfactor(ŷ, Forward,  g, ω);
    s_vector_y_b = create_sfactor(ŷ, Backward, g, ω);

    # Fill the 2D space with layers of appropriate s-factors
    Sx_f_2D = zeros(Complex128, size(g));
    Sx_b_2D = zeros(Complex128, size(g));
    Sy_f_2D = zeros(Complex128, size(g));
    Sy_b_2D = zeros(Complex128, size(g));

    for i = 1:size(g, 1)
        Sy_f_2D[i, :] = s_vector_y_f.^-1;
        Sy_b_2D[i, :] = s_vector_y_b.^-1;
    end

    for i = 1:size(g, 2)
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
