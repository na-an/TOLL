	
    volVectorField adjointTransposeConvection((rho*fvc::grad(U) & Ua));
			
    tmp<fvVectorMatrix> tUaEqn
    (
        rho*fvm::div(-phi, Ua)
      + adjointTransposeConvection
      - rho*fvm::laplacian(nua,Ua)
      + fvm::Sp(alphasolid, Ua)
      + changealpha*adalpha*fvc::grad(alphac)
      ==
        fvOptions(rho,Ua)
    );
	
    fvVectorMatrix& UaEqn = tUaEqn.ref();

    UaEqn.relax();
    fvOptions.constrain(UaEqn);
   
    Uaresidual = solve(UaEqn==-fvc::grad(pa)).initialResidual();
    Info<< "UaEqn.residual: " << Uaresidual << endl;

    Uaresidualsave=Uaresidual[0];

    fvOptions.correct(Ua);
