using Pardiso
using Dates

function dolinearsolve(A::SparseMatrixCSC, b::Array, matrixsym::FDFDMatSymmetry)
    time0 = now()

    @info "Solving linear system"
    @info "Unknowns: $(length(b))"

    if haskey(ENV, "FDFD_SOLVER") && lowercase(ENV["FDFD_SOLVER"]) == "pardiso"
        @info "Solver: Pardiso"
        if haskey(ENV, "OMP_NUM_THREADS")
            @info "Threads: $(ENV["OMP_NUM_THREADS"])"
        else
            @warn "OMP_NUM_THREADS was not set. Pardiso may not be using the expected number of cores."
        end

        # Prefer the MKL Pardiso solver over the vanilla Pardiso solver
        ps = nothing
        try
            ps = MKLPardisoSolver()
        catch e
            if e isa LoadError
                ps = PardisoSolver()
            else
                error("Error while loading the Pardiso solver.")
            end
        end
        set_matrixtype!(ps, Pardiso.COMPLEX_NONSYM)
        x = Pardiso.solve(ps, A, b);
        set_phase!(ps, Pardiso.RELEASE_ALL)
        pardiso(ps)
    else
        @info "Solver: Julia"
        x = lu(A)\b
    end

    @info "Solve time: $(now()-time0)"

    return x
end
