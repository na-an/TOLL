

    volScalarField rAUa(1.0/UaEqn.A());
	volVectorField HbyAa("HbyAa", Ua);
    HbyAa = rAUa*UaEqn.H();
    tUaEqn.clear();

    surfaceScalarField phiHbyAa("phiHbyAa", fvc::flux(HbyAa));

    adjustPhi(phiHbyAa, Ua, pa);

    tmp<volScalarField> rAtUa(rAUa);

    if (simple.consistent())
    {
        rAtUa = 1.0/(1.0/rAUa - UaEqn.H1());
        phiHbyAa +=
            fvc::interpolate(rAtUa() - rAUa)*fvc::snGrad(pa)*mesh.magSf();
        HbyAa -= (rAUa - rAtUa())*fvc::grad(pa);
    }

    tUaEqn.clear();

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(pa, Ua, phiHbyAa, rAtUa());

    // Non-orthogonal pressure corrector loop
    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix paEqn
        (
            fvm::laplacian(rAtUa, pa) == fvc::div(phiHbyAa)
        );

        paEqn.setReference(paRefCell, paRefValue);

        paEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            phia = phiHbyAa - paEqn.flux();
        }
    }


    // Explicitly relax pressure for momentum corrector
    pa.relax();

    // Momentum corrector
    Ua = HbyAa - rAtUa*fvc::grad(pa);
    Ua.correctBoundaryConditions();
    fvOptions.correct(Ua);
    
