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

\*---------------------------------------------------------------------------*/

#include "thermalDiameter.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
    defineTypeNameAndDebug(thermal, 0);

    addToRunTimeSelectionTable
    (
        diameterModel,
        thermal,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::thermal::thermal
(
    const dictionary& diameterProperties,
    const phaseModel& phase
)
:
    diameterModel(diameterProperties, phase),
    dmin_("dmin", dimLength, diameterProperties_),
    dmax_("dmax", dimLength, diameterProperties_),
    Tmin_("Tmin", dimTemperature, diameterProperties_),
    Tmax_("Tmax", dimTemperature, diameterProperties_),
    liquidPhase_(diameterProperties.lookup("liquidPhaseTemperature")),
    d_
    (
        IOobject
        (
            IOobject::groupName("d", phase.name()),
            phase_.time().timeName(),
            phase_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phase_.mesh(),
        dimensionedScalar("d", dimLength, 0.0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::thermal::~thermal()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::diameterModels::thermal::d() const
{
    return d_;
}


void Foam::diameterModels::thermal::correct()
{
    const volScalarField& Tliq = phase_.db().lookupObject<volScalarField>
    (
        liquidPhase_
    );


    forAll (Tliq.internalField(), i)
    {
        if (Tliq[i] < Tmin_.value())
        {
            d_[i] = dmin_.value();
        }
        else if (Tliq[i] > Tmin_.value() &&  Tliq[i] > Tmax_.value())
        {
            d_[i] =
            (
                dmin_.value()*(Tliq[i] - Tmax_.value())
               +dmax_.value()*(Tmin_.value() - Tliq[i])
            )
            /(Tmin_.value()-Tmax_.value());
        }
        else if (Tliq[i] > Tmax_.value())
        {
            d_[i] = dmax_.value();
        }
    }
}


bool Foam::diameterModels::thermal::read(const dictionary& phaseProperties)
{
    diameterModel::read(phaseProperties);

    diameterProperties_.lookup("dmin") >> dmin_;
    diameterProperties_.lookup("dmax") >> dmax_;
    diameterProperties_.lookup("Tmin") >> Tmin_;
    diameterProperties_.lookup("Tmax") >> Tmax_;

    return true;
}


// ************************************************************************* //
