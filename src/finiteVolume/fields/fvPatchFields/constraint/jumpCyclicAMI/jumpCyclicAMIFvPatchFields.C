/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
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

#include "jumpCyclicAMIFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makePatchFieldTypeNames(jumpCyclicAMI);

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
void Foam::jumpCyclicAMIFvPatchField<scalar>::updateInterfaceMatrix
(
    scalarField& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const scalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes
) const
{
    //const labelUList& nbrFaceCells =
    //    this->cyclicAMIPatch().cyclicAMIPatch().neighbPatch().faceCells();

    label nbrPathid = this->cyclicAMIPatch().neighbPatchID();
    const labelUList& nbrFaceCells = lduAddr.patchAddr(nbrPathid);

    scalarField pnf(psiInternal, nbrFaceCells);

    pnf = this->cyclicAMIPatch().interpolate(pnf);

    // only apply jump to original field
    if (&psiInternal == &this->primitiveField())
    {
        Field<scalar> jf(this->jump());

        if (!this->cyclicAMIPatch().owner())
        {
            jf *= -1.0;
        }

        pnf -= jf;
    }

    // Transform according to the transformation tensors
    this->transformCoupleField(pnf, cmpt);

    const labelUList& faceCells = lduAddr.patchAddr(patchId);


    // Multiply the field by coefficients and add into the result
    this->addToInternalField(result, !add, faceCells, coeffs, pnf);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
