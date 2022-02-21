/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    twoLiquidMixingFoam

Description
    Solver for mixing 2 incompressible fluids.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "MULES.H"
#include "subCycle.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"
#include <math.h>
#include <fstream>
#include <iostream>
#include <iosfwd>
#include <stdio.h>
#include "petsc.h"
#include "petscvec.h"
#include <MMA.h>
#include <mpi.h>
template<class Type>
void zeroCells
(
    GeometricField<Type, fvPatchField, volMesh>& vf,
    const labelList& cells
)
{
    forAll(cells, i)
    {
        vf[cells[i]] = Zero;
    }
}

template<class Type>
void oneCells
(
    GeometricField<Type, fvPatchField, volMesh>& vf,
    const labelList& cells
)
{
    forAll(cells, i)
    {
        vf[cells[i]] = 1;
    }
}
template<class Type>
void setCells
(
    GeometricField<Type, fvPatchField, volMesh>& vf,
    const labelList& cells,
    double value
)
{
    forAll(cells, i)
    {
        vf[cells[i]] = value;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

static char help[] = "topology optimization \n";

int main(int argc, char *argv[])
{
    PetscInitialize(&argc,&argv,PETSC_NULL,help);
    
    #include "initContinuityErrs.H"

    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "initContinuityErrs.H"
    
    #include "createFields.H"
    
    #include "initialize_dvg.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    	
    while (simple.loop(runTime))
    {	
        Info<< "Time = " << runTime.timeName() << nl << endl;


        #include "primal_equation.H"

        #include "adjoint_equation.H"

        #include "disspower_equation.H"
        
        #include "sensitivity.H"
        
        #include "costFunction.H"


        if(runTime.writeTime())
        {
            dvg.write();
            dvgPDE.write();
            dvgHea.write();
            
            p.write();
            U.write();
            pa.write();
            Ua.write();
            pc.write();
            Uc.write();
            
            alphac.write();
            adalpha.write();
            
            gsens.write();
            gsensPDE.write();
            gsensHea.write();
            
            fsens.write();
            fsensPDE.write();
            fsensHea.write();

        }
        if ((dvgchange<1e-8) &&(mag(jmixdegree/J1-1.0)<0.01))
        {
            dvg.write();
            dvgPDE.write();
            dvgHea.write();
            
            p.write();
            U.write();
            pa.write();
            Ua.write();
            pc.write();
            Uc.write();
            
            alphac.write();
            adalpha.write();
            
            gsens.write();
            gsensPDE.write();
            gsensHea.write();
            
            fsens.write();
            fsensPDE.write();
            fsensHea.write();
            break;
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
}
    delete mma;
    //delete filter;
    PetscFinalize();
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
