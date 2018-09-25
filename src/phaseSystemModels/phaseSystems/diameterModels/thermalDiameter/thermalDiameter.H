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

Class
    Foam::diameterModels::thermal

Description
    Isothermal dispersed-phase particle diameter model.

SourceFiles
    thermal.C

\*---------------------------------------------------------------------------*/

#ifndef thermalDiameter_H
#define thermalDiameter_H

#include "diameterModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{

/*---------------------------------------------------------------------------*\
                           Class thermal Declaration
\*---------------------------------------------------------------------------*/

class thermal
:
    public diameterModel
{
    // Private data

        //- Reference min diameter for sub-cooling liquid temperature
        dimensionedScalar dmin_;

        //- Reference max diameter
        dimensionedScalar dmax_;

        //- Reference min sub-cooling temperature
        dimensionedScalar Tmin_;

        //- Reference max sub-cooling temperature
        dimensionedScalar Tmax_;

        //- Name of the liquid temperature field
        word liquidPhase_;

        //- Actual diameter field
        volScalarField d_;


public:

    //- Runtime type information
    TypeName("thermal");


    // Constructors

        //- Construct from components
        thermal
        (
            const dictionary& diameterProperties,
            const phaseModel& phase
        );


    //- Destructor
    virtual ~thermal();


    // Member Functions

        //- Return the diameter field
        virtual tmp<volScalarField> d() const;

        //- Correct the diameter field
        virtual void correct();

        //- Read phaseProperties dictionary
        virtual bool read(const dictionary& phaseProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
