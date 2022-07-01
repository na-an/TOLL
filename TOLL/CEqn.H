
    fvScalarMatrix CEqn
    (
        fvm::div(phi,C)
      - fvm::laplacian(Da,C)
      ==
        fvOptions(C)
    );

    CEqn.relax();
    fvOptions.constrain(CEqn);
    CEqn.solve();
    fvOptions.correct(C);

Info<< "C volume fraction = "
        << C.weightedAverage(mesh.Vsc()).value()
        << "  Min(C) = " << min(C).value()
        << "  Max(C) = " << max(C).value()
        << endl;