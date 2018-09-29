/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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

#include "foamVtkFieldTypes.H"

#include "labelIOField.H"
#include "scalarIOField.H"
#include "sphericalTensorIOField.H"
#include "symmTensorIOField.H"
#include "tensorIOField.H"

#include "areaFields.H"
#include "volFields.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::wordList Foam::vtk::fieldTypes::vols
({
    Foam::volScalarField::typeName,
    Foam::volVectorField::typeName,
    Foam::volSphericalTensorField::typeName,
    Foam::volSymmTensorField::typeName,
    Foam::volTensorField::typeName
});


const Foam::wordList Foam::vtk::fieldTypes::dims
({
    Foam::volScalarField::Internal::typeName,
    Foam::volVectorField::Internal::typeName,
    Foam::volSphericalTensorField::Internal::typeName,
    Foam::volSymmTensorField::Internal::typeName,
    Foam::volTensorField::Internal::typeName
});


const Foam::wordList Foam::vtk::fieldTypes::points
({
    Foam::pointScalarField::typeName,
    Foam::pointVectorField::typeName,
    Foam::pointSphericalTensorField::typeName,
    Foam::pointSymmTensorField::typeName,
    Foam::pointTensorField::typeName
});


const Foam::wordList Foam::vtk::fieldTypes::clouds
({
    Foam::labelIOField::typeName,
    Foam::scalarIOField::typeName,
    Foam::vectorIOField::typeName,
    Foam::symmTensorIOField::typeName,
    Foam::tensorIOField::typeName
});


const Foam::wordList Foam::vtk::fieldTypes::areas
({
    Foam::areaScalarField::typeName,
    Foam::areaVectorField::typeName,
    Foam::areaSphericalTensorField::typeName,
    Foam::areaSymmTensorField::typeName,
    Foam::areaTensorField::typeName
});


// ************************************************************************* //
