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

#include "CaOutletFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "turbulentTransportModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CaOutletFvPatchScalarField::
CaOutletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{}


Foam::CaOutletFvPatchScalarField::
CaOutletFvPatchScalarField
(
    const CaOutletFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::CaOutletFvPatchScalarField::
CaOutletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict)
{}


Foam::CaOutletFvPatchScalarField::
CaOutletFvPatchScalarField
(
    const CaOutletFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::CaOutletFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchField<scalar>& Cp =
        patch().lookupPatchField<volScalarField, scalar>("C");
	const fvPatchField<scalar>& Cap =
        patch().lookupPatchField<volScalarField, scalar>("Ca");		
    const fvPatchField<scalar>& Cvp =
        patch().lookupPatchField<volScalarField, scalar>("Cv");
    const fvPatchField<vector>& Up =
        patch().lookupPatchField<volVectorField, vector>("U");
    const dictionary& transportProperties = db().lookupObject<IOdictionary>("transportProperties");
    dimensionedScalar Dap(transportProperties.lookup("Da"));
    const scalarField& deltainv = patch().deltaCoeffs();

    scalarField Up_n = (Up & patch().nf());
    scalarField Caneigh_n = Cap.patchInternalField();

    operator == ((Dap.value()*Caneigh_n*deltainv-2*(Cp-Cvp))/(Dap.value()*deltainv+Up_n));

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::CaOutletFvPatchScalarField::write(Ostream& os) const
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
        CaOutletFvPatchScalarField
    );
}

// ************************************************************************* //
