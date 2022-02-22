
    volScalarField rAUc(1.0/UcEqn.A());
    volVectorField HbyAc("HbyAc", Uc);
    HbyAc = rAUc*UcEqn.H();
    tUcEqn.clear();
    surfaceScalarField phiHbyAc("phiHbyAc", fvc::flux(HbyAc));
	
    adjustPhi(phiHbyAc, Uc, pc);

    tmp<volScalarField> rAtUc(rAUc);

    if (simple.consistent())
    {
        rAtUc = 1.0/(1.0/rAUc - UcEqn.H1());
        phiHbyAc +=
            fvc::interpolate(rAtUc() - rAUc)*fvc::snGrad(pc)*mesh.magSf();
        HbyAc -= (rAUc - rAtUc())*fvc::grad(pc);
    }

    tUcEqn.clear();
    // Update the pressure BCs to ensure flux consistency
    constrainPressure(pc, Uc, phiHbyAc, rAtUc());

    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pcEqn
        (
            fvm::laplacian(rAtUc, pc) == fvc::div(phiHbyAc)
        );

        pcEqn.setReference(pcRefCell, pcRefValue);

        pcEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            phic = phiHbyAc - pcEqn.flux();


        }
    }

    // Explicitly relax pressure for momentum corrector
    pc.relax();

    // Momentum corrector
    Uc = HbyAc - rAtUc()*fvc::grad(pc);
    Uc.correctBoundaryConditions();
    fvOptions.correct(Uc);