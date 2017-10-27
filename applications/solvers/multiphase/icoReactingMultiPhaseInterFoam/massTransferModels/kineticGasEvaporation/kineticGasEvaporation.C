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

#include "kineticGasEvaporation.H"
#include "constants.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcDiv.H"
#include "surfaceInterpolate.H"
#include "fvcReconstruct.H"
#include "fvm.H"
#include "zeroGradientFvPatchFields.H"

using namespace Foam::constant;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::meltingEvaporationModels::kineticGasEvaporation<Thermo, OtherThermo>
::kineticGasEvaporation
(
    const dictionary& dict,
    const phasePair& pair
)
:
    InterfaceCompositionModel<Thermo, OtherThermo>(dict, pair),
    C_("C",  dimless, dict.lookup("C")),
    Tactivate_
    (
        "Tactivate",
        dimTemperature,
        dict.lookup("Tactivate")
    ),
    Mv_
    (
        "Mv",
        dimMass/dimMoles,
        dict.lookupOrDefault<scalar>("Mv", -1)
    ),
    phi0_
    (
        IOobject
        (
            "phi0",
            this->mesh_.time().timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("zero", dimMass/dimTime/dimVolume, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    alphaMax_(dict.lookupOrDefault<scalar>("alphaMax", 0.03)),
    alphaMin_(dict.lookupOrDefault<scalar>("alphaMin", 0.01))
{
    word fullSpeciesName = this->transferSpecie();

    const label tempOpen(fullSpeciesName.find('.'));

    const word speciesName(fullSpeciesName(0, tempOpen));

    // Get the continuous thermo
    const typename OtherThermo::thermoType& localThermo =
        this->getLocalThermo
        (
            speciesName,
            this->toThermo_
        );

    Mv_.value() = localThermo.W();

    if (Mv_.value() == -1)
    {
        FatalErrorIn
        (
           "meltingEvaporationModels::"
           "kineticGasEvaporation<Thermo, OtherThermo>::"
           "kineticGasEvaporation"
           "("
            "const dictionary& dict,"
            "const phasePair& pair"
           ")"
        )
            << " Please provide the molar weight (Mv) of vapour [g/mol] "
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::meltingEvaporationModels::kineticGasEvaporation<Thermo, OtherThermo>
::Kexp(label variable, const volScalarField& field)
{

    if (this->modelVariable_ == variable)
    {
        const volScalarField& to = this->pair().to();

        const volScalarField& from = this->pair().from();

        const fvMesh& mesh = this->mesh_;

        const volScalarField& T = mesh.lookupObject<volScalarField>("T").oldTime();

        const dimensionedScalar HerztKnudsConst
        (
            sqrt
            (
                Mv_
               /2.0
               /constant::physicoChemical::R
               /mathematical::pi
               /pow3(Tactivate_)
            )
        );

        word fullSpeciesName = this->transferSpecie();

        const label tempOpen(fullSpeciesName.find('.'));

        const word speciesName(fullSpeciesName(0, tempOpen));

        tmp<volScalarField> L = this->L(speciesName, field);

        const volVectorField gradFrom(fvc::grad(from));
        const volVectorField gradTo(fvc::grad(to));

        volScalarField areaDensity
        (
            "areaDensity",
            sqrt(mag(gradFrom)*mag(gradTo))
        );


        volScalarField gradAlphaf
        (
            gradFrom
          & gradTo
        );

//         const volScalarField trigger
//         (
//             pos(alphaMax_ - from)
//         );

       // const volScalarField mask(nearInterface*gradAlphaf);

//         const volScalarField Tmask
//         (
//             "Tmask",
//             pos(alphaMax_ - from)*pos(from - alphaMin_)
//         );

        volScalarField Tmask
        (
            "Tmask",
            neg(from + 1.0)
        );

        Tmask *= 0.0;

        forAll (Tmask, celli)
        {
            if (gradAlphaf[celli] < 0)
            {
                if (from[celli] > alphaMin_ &&  from[celli] < alphaMax_)
                {
                    scalar alphaRes = 1.0 - from[celli] - to[celli];
                    if (alphaRes < 0.01)
                    {
                        Tmask[celli] = 1.0;
                    }
                }
            }
        }

        const volScalarField Tave(T*Tmask);

        volScalarField massFluxEvap
        (
            "massFluxEvap",
            2*C_/(2 - C_)
          * HerztKnudsConst
          * L()
          * this->pair().to().rho()
          * max
            (
                (Tave - Tactivate_),
                dimensionedScalar("T0", dimTemperature, 0.0)
            )
        );

         // Liquid normalization
        const dimensionedScalar Nl
        (
            gSum((areaDensity*mesh.V())())
           /(
               gSum
               (
                   ((areaDensity*from)*mesh.V())()
               )
             + dimensionedScalar("SMALL", inv(dimless), VSMALL)
            )
        );

        // Local density rate (kg/m3/s)
        phi0_ = massFluxEvap*Nl*areaDensity*from; //

        if (this->pair().from().mesh().time().outputTime())
        {
            areaDensity.write();
            Tmask.write();
        }

        return (1.0*phi0_);
    }
    else
    {
        return tmp<volScalarField> ();
    }
}


template<class Thermo, class OtherThermo>
const Foam::dimensionedScalar&
Foam::meltingEvaporationModels::kineticGasEvaporation<Thermo, OtherThermo>
::Tactivate() const
{
    return Tactivate_;
}


// ************************************************************************* //
