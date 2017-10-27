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

Class
    Foam::phaseChangeTwoPhaseMixtures::kineticGasEvaporation

Description
    kineticGasEvaporation evaporation model.

SourceFiles
    kineticGasEvaporation.C

\*---------------------------------------------------------------------------*/

#ifndef kineticGasEvaporation_H
#define kineticGasEvaporation_H


#include "InterfaceCompositionModel.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * * * *//

namespace Foam
{

class phasePair;

namespace meltingEvaporationModels
{

/*---------------------------------------------------------------------------*\
                   Class kineticGasEvaporation
\*---------------------------------------------------------------------------*/

template<class Thermo, class OtherThermo>
class kineticGasEvaporation
:
    public InterfaceCompositionModel<Thermo, OtherThermo>
{
    // Private data

        //- Evaporation coefficient
        dimensionedScalar C_;

        //- Activation temperature
        const dimensionedScalar Tactivate_;

        //- Molar weight of the vapour in the continous phase
        dimensionedScalar Mv_;

        //- Source
        volScalarField phi0_;

        //- Alpha maximum value for the mass transfer
        scalar alphaMax_;

        //- Alpha minumum value for the mass transfer
        scalar alphaMin_;


public:

    //- Runtime type information
    TypeName("kineticGasEvaporation");


    // Constructors

        //- Construct from components
        kineticGasEvaporation
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~kineticGasEvaporation()
    {}


    // Member Functions

        //- Explicit mass transfer coefficient
        virtual tmp<volScalarField> Kexp
        (
            label variable,
            const volScalarField& field
        );

        //- Return Tref
        virtual const dimensionedScalar& Tactivate() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace meltingEvaporationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "kineticGasEvaporation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
