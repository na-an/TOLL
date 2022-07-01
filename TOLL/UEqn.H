
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
        fvOptions.correct(U);
    }
    Uresidual0 = Uresidual[0];
