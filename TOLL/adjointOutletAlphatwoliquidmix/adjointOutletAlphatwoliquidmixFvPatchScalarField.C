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

#include "adjointOutletAlphatwoliquidmixFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "turbulentTransportModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adjointOutletAlphatwoliquidmixFvPatchScalarField::
adjointOutletAlphatwoliquidmixFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{}


Foam::adjointOutletAlphatwoliquidmixFvPatchScalarField::
adjointOutletAlphatwoliquidmixFvPatchScalarField
(
    const adjointOutletAlphatwoliquidmixFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::adjointOutletAlphatwoliquidmixFvPatchScalarField::
adjointOutletAlphatwoliquidmixFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict)
{}


Foam::adjointOutletAlphatwoliquidmixFvPatchScalarField::
adjointOutletAlphatwoliquidmixFvPatchScalarField
(
    const adjointOutletAlphatwoliquidmixFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::adjointOutletAlphatwoliquidmixFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

	const fvPatchField<scalar>& alphacp =
        patch().lookupPatchField<volScalarField, scalar>("alphac");
	const fvPatchField<scalar>& adalphap =
        patch().lookupPatchField<volScalarField, scalar>("adalpha");		
    const fvPatchField<scalar>& alphadp =
        patch().lookupPatchField<volScalarField, scalar>("alphad");
    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>("phi");

    const fvsPatchField<scalar>& phiap =
        patch().lookupPatchField<surfaceScalarField, scalar>("phia");

    const fvPatchField<vector>& Up =
        patch().lookupPatchField<volVectorField, vector>("U");

    const fvPatchField<vector>& Uap =
        patch().lookupPatchField<volVectorField, vector>("Ua");

    //scalarField Up_n = phip / patch().magSf();
    scalarField Up_n = (Up & patch().nf());
    scalarField Uap_n = phiap / patch().magSf();

	const dictionary& transportProperties = db().lookupObject<IOdictionary>("transportProperties");
    dimensionedScalar Dabp(transportProperties.lookup("Dabbound"));

    dimensionedScalar nua(transportProperties.lookup("nua"));

    const scalarField& deltainv = patch().deltaCoeffs();
	
    scalarField adalphaneigh_n = adalphap.patchInternalField();
	
    operator == ((Dabp.value()*adalphaneigh_n*deltainv-2*(alphacp-alphadp))/(Dabp.value()*deltainv+Up_n));

	fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::adjointOutletAlphatwoliquidmixFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        adjointOutletAlphatwoliquidmixFvPatchScalarField
    );
}

// ************************************************************************* //
