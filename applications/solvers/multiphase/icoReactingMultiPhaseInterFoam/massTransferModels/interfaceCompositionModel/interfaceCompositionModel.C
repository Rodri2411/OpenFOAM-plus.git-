/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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

#include "interfaceCompositionModel.H"
#include "phaseModel.H"
#include "phasePair.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(interfaceCompositionModel, 0);
    defineRunTimeSelectionTable(interfaceCompositionModel, dictionary);
}


namespace Foam
{
    template<>
    const char* Foam::NamedEnum
    <
        Foam::interfaceCompositionModel::modelVariable,
        3
    >::names[] =
    {
        "temperature",
        "pressure",
        "massFraction"
    };


} // End namespace Foam

const Foam::NamedEnum
<
    Foam::interfaceCompositionModel::modelVariable,
    3
> Foam::interfaceCompositionModel::modelVariableNames;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceCompositionModel::interfaceCompositionModel
(
    const dictionary& dict,
    const phasePair& pair
)
:
    modelVariable_(modelVariableNames.read(dict.lookup("variable"))),
    semiImplicit_(dict.lookupOrDefault<bool>("semiImplicit", false)),
    pair_(pair),
    speciesName_(dict.lookupOrDefault<word>("species", "none")),
    mesh_(pair_.from().mesh())
{/*
    if (dict.found("species"))
    {
        speciesNames_ = ;
    }*/
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interfaceCompositionModel::~interfaceCompositionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::word Foam::interfaceCompositionModel::transferSpecie() const
{
    return speciesName_;
}


const Foam::phasePair& Foam::interfaceCompositionModel::pair() const
{
    return pair_;
}


const Foam::word Foam::interfaceCompositionModel::variable() const
{
    return modelVariableNames[modelVariable_];
}


bool Foam::interfaceCompositionModel::semiImplicit() const
{
    return semiImplicit_;
}


Foam::tmp<Foam::volScalarField> Foam::interfaceCompositionModel::KexpEnergy
(
    label modelVariable,
    const volScalarField& field
)
const
{
    NotImplemented;
    return tmp<volScalarField>();
}


// ************************************************************************* //
