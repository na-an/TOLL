
    tmp<fvVectorMatrix> tUEqn
    (

        rho*fvm::div(phi, U)
      + rho*turbulence->divDevReff(U)
      + fvm::Sp(alphasolid, U)
      ==
        fvOptions(rho,U)
    );

    
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (simple.momentumPredictor())
    {
        Uresidual = solve(UEqn == -fvc::grad(p)).initialResidual();
        Info<< "UEqn.residual: " << Uresidual << endl;

        Uresidualsave=Uresidual[0];
        fvOptions.correct(U);

    }
    
Info<< "UEqn.residual: " << Uresidualsave << endl;

