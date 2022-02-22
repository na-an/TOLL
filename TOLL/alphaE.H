
    fvScalarMatrix alphacE
    (
       fvm::div(phi,alphac)
     - fvm::laplacian(Dab,alphac)
	 ==
       fvOptions(alphac)
    );


    alphacE.relax();
    fvOptions.constrain(alphacE);
    alphacE.solve();
    fvOptions.correct(alphac);


Info<< "alpha volume fraction = "
        << alphac.weightedAverage(mesh.Vsc()).value()
        << "  Min() = " << min(alphac).value()
        << "  Max() = " << max(alphac).value()
        << endl;