using Pardiso, MUMPS

function dolinearsolve(A::SparseMatrixCSC, b::Array, matrixsym::FDFDMatSymmetry)
    tic();

    print_info("Solving linear system")
    print_info(@sprintf("Num unknowns: %.2E", length(b)))

    haskey(ENV, "OMP_NUM_THREADS") && print_info("Num threads: ", ENV["OMP_NUM_THREADS"])

    if haskey(ENV, "FDFD_SOLVER") && lowercase(ENV["FDFD_SOLVER"]) == "pardiso"
        print_info("Solver: Pardiso")
        ps = PardisoSolver()

        matrixsym == CNSym ? sym = Pardiso.COMPLEX_NONSYM : sym = Pardiso.COMPLEX_SYM
        set_matrixtype!(ps, sym)

        A_pardiso = get_matrix(ps, A, :N)

        set_phase!(ps, Pardiso.ANALYSIS)
        pardiso(ps, A_pardiso, b)

        set_phase!(ps, Pardiso.NUM_FACT)
        pardiso(ps, A_pardiso, b)

        set_phase!(ps, Pardiso.SOLVE_ITERATIVE_REFINE)
        x = similar(b)
        pardiso(ps, x, A_pardiso, b)

        set_phase!(ps, Pardiso.RELEASE_ALL)
        pardiso(ps)
    # elseif haskey(ENV, "FDFD_SOLVER") && lowercase(ENV["FDFD_SOLVER"]) == "mumps"
    #     print_info("Solver: MUMPS")
    #     x = similar(b)
    #
    #     matrixsym == CNSym ? sym = 0 : sym = 2
    #     Afact = factorMUMPS(A, sym)
    #
    #     x = applyMUMPS(Afact, b)
    #
    #     destroyMUMPS(Afact)
    else
        print_info("Solver: Julia")
        x = lufact(A)\b
    end

    print_info(@sprintf("Solve time: %.2f minutes", toq()/60))

    return x
end
