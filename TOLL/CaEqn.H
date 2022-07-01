    fvScalarMatrix CaEqn
    (
      - fvm::div(-phi,Ca)
      + fvm::laplacian(Da,Ca)
      ==
        fvOptions(Ca)
    );

    CaEqn.relax();
    fvOptions.constrain(CaEqn);
    CaEqn.solve();
    fvOptions.correct(Ca);

Info<< "Ca average value = "
        << Ca.weightedAverage(mesh.Vsc()).value()
        << "  Min(Ca) = " << min(Ca).value()
        << "  Max(Ca) = " << max(Ca).value()
        << endl;