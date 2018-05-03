function dolinearsolve(A::SparseMatrixCSC, b::Array; matrixtype=Pardiso.COMPLEX_NONSYM)
    tic();
    print_info("Solving linear system");
    print_info(@sprintf("Num unknowns: %.2E", length(b)));
    pardiso_success = false;
    try
        ps = PardisoSolver();
        set_matrixtype!(ps, matrixtype);
        set_solver!(ps, Pardiso.DIRECT_SOLVER);
        pardisoinit(ps);
        try print_info("Num threads: ", ENV["OMP_NUM_THREADS"]) end
        print_info("Solver: Pardiso");
        x = Pardiso.solve(ps, A, b);
        pardiso_success = true;
    end
    if ~pardiso_success
        print_info("Solver: Julia");
        x = lufact(A)\b;
    end
    print_info(@sprintf("Solve time: %.2f minutes", toq()/60));
    return x
end