/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "adjointOutletVelocitytwoliquidmixFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvPatchFieldMapper.H"
#include "turbulentTransportModel.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adjointOutletVelocitytwoliquidmixFvPatchVectorField::
adjointOutletVelocitytwoliquidmixFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF)
{}


Foam::adjointOutletVelocitytwoliquidmixFvPatchVectorField::
adjointOutletVelocitytwoliquidmixFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict)
{}


Foam::adjointOutletVelocitytwoliquidmixFvPatchVectorField::
adjointOutletVelocitytwoliquidmixFvPatchVectorField
(
    const adjointOutletVelocitytwoliquidmixFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper)
{}


Foam::adjointOutletVelocitytwoliquidmixFvPatchVectorField::
adjointOutletVelocitytwoliquidmixFvPatchVectorField
(
    const adjointOutletVelocitytwoliquidmixFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update the coefficients associated with the patch field
void Foam::adjointOutletVelocitytwoliquidmixFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    const fvPatchField<vector>& Up =
        patch().lookupPatchField<volVectorField, vector>("U");

    const fvPatchField<vector>& Uap =
        patch().lookupPatchField<volVectorField, vector>("Ua");
		
    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>("phi");   
		
	const fvsPatchField<scalar>& phiap =
        patch().lookupPatchField<surfaceScalarField, scalar>("phia");

    const dictionary& transportProperties = db().lookupObject<IOdictionary>("transportProperties");
    dimensionedScalar nua(transportProperties.lookup("nua"));

    const scalarField& deltainv = patch().deltaCoeffs();  // dist^(-1)

    //Primal velocity , mag of normal component and tangential component 

    scalarField Up_ns = phip / patch().magSf(); //Mag. of normal

    vectorField Up_t = Up - (phip * patch().Sf())/(patch().magSf()*patch().magSf()); //Tangential


    // Neighbouring node's vel.
    vectorField Uaneigh = Uap.patchInternalField();
    vectorField Uaneigh_n = (Uaneigh & patch().nf())*patch().nf(); //Normal
    vectorField Uaneigh_t = Uaneigh - Uaneigh_n; //Tangential   
		
	vectorField Uap_t = (nua.value()*deltainv*Uaneigh_t) / (Up_ns+nua.value()*deltainv) ;

    vectorField Uap_n = (phiap * patch().Sf())/(patch().magSf()*patch().magSf());
		
    vectorField::operator==(Uap_t + Uap_n);

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::adjointOutletVelocitytwoliquidmixFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        adjointOutletVelocitytwoliquidmixFvPatchVectorField
    );
}


// ************************************************************************* //
