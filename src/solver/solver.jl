using Pardiso

function dolinearsolve(A::SparseMatrixCSC, b::Array, matrixsym::FDFDMatSymmetry)
    tic();

    @info "Solving linear system"
    @info "Num unknowns: $(length(b))"

    haskey(ENV, "OMP_NUM_THREADS") && @info "Num threads: " ENV["OMP_NUM_THREADS"]

    if haskey(ENV, "FDFD_SOLVER") && lowercase(ENV["FDFD_SOLVER"]) == "pardiso"
        @info "Solver: Pardiso"

        ps = PardisoSolver()
        matrixsym == CNSym ? sym = Pardiso.COMPLEX_NONSYM : sym = Pardiso.COMPLEX_SYM
        set_matrixtype!(ps, sym)
        x = Pardiso.solve(ps, A, b);
        set_phase!(ps, Pardiso.RELEASE_ALL)
        pardiso(ps)
    else
        @info "Solver: Julia"

        x = lufact(A)\b
    end

    @info "Solve time: $(toq()/60) minutes"

    return x
end
