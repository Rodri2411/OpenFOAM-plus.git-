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
    Foam::phaseChangeTwoPhaseMixtures::Lee

Description
    Mass tranfer Lee model. Simple model driven by field value difference as:

    mDot = C*rho*alpha*(value - Tactivate)/Tactivate

    if C > 0:
        mDot = C*rho*alpha*(value - Tactivate)/Tactivate for value > Tactivate
        mDot = 0.0                                           value < Tactivate

    if C < 0:
        mDot = -C*rho*alpha*(Tactivate - value)/Tactivate for value < Tactivate
        mDot = 0.0                                            value > Tactivate


SourceFiles
    Lee.C

\*---------------------------------------------------------------------------*/

#ifndef Lee_H
#define Lee_H

#include "InterfaceCompositionModel.H"


// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * * * *//

namespace Foam
{
namespace meltingEvaporationModels
{

/*---------------------------------------------------------------------------*\
                           Class Lee
\*---------------------------------------------------------------------------*/

template<class Thermo, class OtherThermo>
class Lee
:
    public InterfaceCompositionModel<Thermo, OtherThermo>
{
    // Private data

        //- Condensation coefficient [1/s]
        dimensionedScalar C_;

        //- Phase transition temperature
        const dimensionedScalar Tactivate_;

        //- Phase minimum value for activation
        scalar alphaMin_;


public:

    //- Runtime type information
    TypeName("Lee");


    // Constructors

        //- Construct from components
        Lee
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~Lee()
    {}


    // Member Functions

        //- Explicit mass transfer coefficient
        virtual tmp<volScalarField> Kexp
        (
            label variable,
            const volScalarField& field
        ) const;

        //- Semi implicit species mass transfer coefficient
        virtual tmp<volScalarField> Kimp
        (
            label variable,
            const volScalarField& field
        ) const;

        //- Return T transition between phases
        virtual const dimensionedScalar& Tactivate() const;

        //- Derivate sign of the source term:
        // d(Sk-i)/d(variable) > 0 or d(Si-k)/d(variable) < 0
        virtual label dSdVariable();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace meltingEvaporationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "Lee.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
