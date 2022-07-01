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

#include "paOutletFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "turbulentTransportModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::paOutletFvPatchScalarField::
paOutletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{}


Foam::paOutletFvPatchScalarField::
paOutletFvPatchScalarField
(
    const paOutletFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::paOutletFvPatchScalarField::
paOutletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict)
{}


Foam::paOutletFvPatchScalarField::
paOutletFvPatchScalarField
(
    const paOutletFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::paOutletFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
		
    const fvPatchField<vector>& Up =
        patch().lookupPatchField<volVectorField, vector>("U");
    const fvPatchField<vector>& Uap =
        patch().lookupPatchField<volVectorField, vector>("Ua");
    const dictionary& transportProperties = db().lookupObject<IOdictionary>("transportProperties");
    dimensionedScalar nup(transportProperties.lookup("nu"));
    const scalarField& deltainv = patch().deltaCoeffs();

    scalarField Up_n = Up & patch().nf();
    scalarField Uap_n = Uap & patch().nf();
    scalarField Uaneigh_n = (Uap.patchInternalField() & patch().nf());

	operator ==((Up_n * Uap_n) + nup.value()*deltainv*(Uap_n-Uaneigh_n));
    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::paOutletFvPatchScalarField::write(Ostream& os) const
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
        paOutletFvPatchScalarField
    );
}

// ************************************************************************* //
