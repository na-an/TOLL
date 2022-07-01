	
    volVectorField adjointTransposeConvection((rho*fvc::grad(U) & Ua));
			
    tmp<fvVectorMatrix> tUaEqn
    (
        rho*fvm::div(-phi, Ua)
      + adjointTransposeConvection
      - rho*fvm::laplacian(nu,Ua)
      + fvm::Sp(alphasolid, Ua)
      + changeC*Ca*fvc::grad(C)
      ==
        fvOptions(rho,Ua)
    );
	
    fvVectorMatrix& UaEqn = tUaEqn.ref();

    UaEqn.relax();
    fvOptions.constrain(UaEqn);
   
    if (simple.momentumPredictor())
    {
        Uaresidual = solve(UaEqn==-fvc::grad(pa)).initialResidual();
        fvOptions.correct(Ua);
    }
    Uaresidual0 = Uaresidual[0];

    fvOptions.correct(Ua);
