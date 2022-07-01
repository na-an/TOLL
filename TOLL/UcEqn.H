			
    volVectorField adjointTransposeConvectionc((rho*fvc::grad(U) & Uc));

    tmp<fvVectorMatrix> tUcEqn
    (
        rho*fvm::div(-phi, Uc)
      + adjointTransposeConvectionc
      - rho*fvm::laplacian(nu,Uc)
      + fvm::Sp(alphasolid, Uc)
      ==
        fvOptions(rho,Uc)
    );
	
    fvVectorMatrix& UcEqn = tUcEqn.ref();

    UcEqn.relax();
    
    fvOptions.constrain(UcEqn);

    if (simple.momentumPredictor())
    {
        Ucresidual = solve(UcEqn == -fvc::grad(pc)).initialResidual();        
        fvOptions.correct(Uc);
    }

    Ucresidual0 = Ucresidual[0];
