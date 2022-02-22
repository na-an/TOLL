
    fvScalarMatrix adalphaE
    (
      - fvm::div(-phi,adalpha)
      + fvm::laplacian(Dab,adalpha)
      ==
        fvOptions(adalpha)
    );

    adalphaE.relax();
    fvOptions.constrain(adalphaE);
    adalphaE.solve();
    fvOptions.correct(adalpha);


Info<< "adalpha volume fraction = "
        << adalpha.weightedAverage(mesh.Vsc()).value()
        << "  Min() = " << min(adalpha).value()
        << "  Max() = " << max(adalpha).value()
        << endl;