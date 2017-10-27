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

#include "MassTransferPhaseSystem.H"

#include "HashPtrTable.H"

#include "fvcDiv.H"
#include "fvmSup.H"
#include "fvMatrix.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::MassTransferPhaseSystem<BasePhaseSystem>::
MassTransferPhaseSystem(const fvMesh& mesh)
:
    BasePhaseSystem(mesh)
{
    this->generatePairsAndSubModels("massTransferModel", massTransferModels_);

    forAllConstIter(massTransferModelTable, massTransferModels_, iterModel)
    {
        if (!dmdt_.found(iterModel()->pair()))
        {
            dmdt_.insert
            (
                iterModel()->pair(),
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName("dmdt", iterModel()->pair().name()),
                        this->mesh().time().timeName(),
                        this->mesh(),
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    this->mesh(),
                    dimensionedScalar("zero", dimDensity/dimTime, 0)
                )
            );

            // Create dmdtYi (explicit mass transfer) for each species
            word specie = iterModel()->transferSpecie();
            word modelVariable = iterModel()->variable();

            // Create variable as: variable + specieName (i.e. T O2)
            // to denote O2 mass transfer driven by T
            Pair<word> variableSpecie(modelVariable, specie);

            dmdtYiTable dmdtYi;

            // Insert volScalarField named : pair() + T + O2 into table
            dmdtYi.insert
            (
                variableSpecie,
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName
                        (
                            "dmdt",
                            iterModel()->pair().name()
                            + variableSpecie.first()
                            + variableSpecie.second()
                        ),
                        this->mesh().time().timeName(),
                        this->mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    this->mesh(),
                    dimensionedScalar("zero", inv(dimTime), 0)
                )
            );

            // Fill mass tranfer per pair + driven variable (p, T or y) + specie
            // name
            dmdtYi_.insert(iterModel()->pair(), dmdtYi);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::MassTransferPhaseSystem<BasePhaseSystem>::
~MassTransferPhaseSystem()
{}

// * * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * //

template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::MassTransferPhaseSystem<BasePhaseSystem>::calculateL
(
    const volScalarField& dmdtNetki,
    const phasePairKey& keyik,
    const phasePairKey& keyki,
    const volScalarField& T
) const
{
    tmp<volScalarField> tL
    (
        new volScalarField
        (
            IOobject
            (
                "tL",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero", dimEnergy/dimMass, 0)
        )
    );
    volScalarField& L = tL.ref();

    if (massTransferModels_.found(keyik))
    {
        const autoPtr<interfaceCompositionModel>& interfacePtr =
            massTransferModels_[keyik];

        word speciesName = interfacePtr->transferSpecie();

        const label tempOpen(speciesName.find('.'));

        const word species(speciesName(0, tempOpen));

        L -= neg(dmdtNetki)*interfacePtr->L(species, T);
    }

    if (massTransferModels_.found(keyki))
    {
        const autoPtr<interfaceCompositionModel>& interfacePtr =
            massTransferModels_[keyki];

        word speciesName = interfacePtr->transferSpecie();

        const label tempOpen(speciesName.find('.'));

        const word species(speciesName(0, tempOpen));

        L += pos(dmdtNetki)*interfacePtr->L(species, T);
    }

    return tL;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::MassTransferPhaseSystem<BasePhaseSystem>::dmdt
(
    const phasePairKey& key
) const
{
    tmp<volScalarField> tdmdt
    (
        new volScalarField
        (
            IOobject
            (
                "dmdt",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero", dimDensity/dimTime, 0)
        )
    );

    volScalarField& dmdt = tdmdt.ref();

    if (dmdt_.found(key))
    {
        dmdt = *dmdt_[key];
    }

    return tdmdt;
}


template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::MassTransferPhaseSystem<BasePhaseSystem>::dmdtYi
(
    const word& phaseSpeciesName
) const
{
    const label tempOpen = phaseSpeciesName.find('.');
    const word species = phaseSpeciesName(0, tempOpen);
    const word phasei = phaseSpeciesName
    (
        tempOpen + 1, phaseSpeciesName.size()
    );

    tmp<volScalarField> tdmdt
    (
        new volScalarField
        (
            IOobject
            (
                "dmdt",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero", dimDensity/dimTime, 0)
        )
    );
    volScalarField& dmdt = tdmdt.ref();

    // Look in all other phases for mass transfer models
    forAllConstIter(phaseSystem::phaseModelTable, this->phaseModels_, iterk)
    {
        const phaseModel& phasek = iterk()();

        if (phasei != iterk()().name())
        {
            // Phase i to phase k
            const phasePairKey keyik
            (
                phasei,
                phasek.name(),
                true
            );

            // Phase k to phase i
            const phasePairKey keyki
            (
                phasek.name(),
                phasei,
                true
            );

            // phase i to Phase k (negative)
            if (dmdtYi_.found(keyik))
            {
                // T based
                Pair<word> TSpecie
                (
                    interfaceCompositionModel::modelVariableNames
                    [
                        interfaceCompositionModel::T
                    ],
                    species
                );

                if (dmdtYi_[keyik].found(TSpecie))
                {
                    dmdt -= *dmdtYi_[keyik][TSpecie];
                }

                // Pressure based
                Pair<word> PSpecie
                (
                    interfaceCompositionModel::modelVariableNames
                    [
                        interfaceCompositionModel::P
                    ],
                    species
                );

                if (dmdtYi_[keyik].found(PSpecie))
                {
                    dmdt -= *dmdtYi_[keyik][PSpecie];
                }
            }

            // phase k to Phase i (positive)
            if (dmdtYi_.found(keyki))
            {
                Pair<word> TSpecie
                (
                    interfaceCompositionModel::modelVariableNames
                    [
                        interfaceCompositionModel::T
                    ],
                    species
                );

                if (dmdtYi_[keyki].found(TSpecie))
                {
                    dmdt += *dmdtYi_[keyki][TSpecie];
                }

                Pair<word> PSpecie
                (
                    interfaceCompositionModel::modelVariableNames
                    [
                        interfaceCompositionModel::P
                    ],
                    species
                );

                if (dmdtYi_[keyki].found(PSpecie))
                {
                    dmdt += *dmdtYi_[keyki][PSpecie];
                }
            }
        }
    }

    return tdmdt;
}


template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::MassTransferPhaseSystem<BasePhaseSystem>::dmdtYi
(
    const phasePairKey& key,
    word var,
    word specieName
) const
{

    Pair<word> variableSpecie(var, specieName);

    tmp<volScalarField> tdmdt
    (
        new volScalarField
        (
            IOobject
            (
                "dmdt",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero", dimDensity/dimTime, 0)
        )
    );

    if (dmdtYi_.found(key))
    {
        return *dmdtYi_[key].operator[](variableSpecie);
    }

    return tdmdt;
}


template<class BasePhaseSystem>
Foam::volScalarField&
Foam::MassTransferPhaseSystem<BasePhaseSystem>::dmdtYi
(
    const phasePairKey& key,
    word var,
    word specieName
)
{
    Pair<word> variableSpecie(var, specieName);
    return *dmdtYi_[key].operator[](variableSpecie);
}


template<class BasePhaseSystem>
Foam::tmp<Foam::fvScalarMatrix>
Foam::MassTransferPhaseSystem<BasePhaseSystem>::heatTransfer
(
    const volScalarField& T
)
{
    tmp<fvScalarMatrix> tEqnPtr
    (
        new fvScalarMatrix(T, dimEnergy/dimTime)
    );

    fvScalarMatrix& eqn = tEqnPtr.ref();

    forAllIter(phaseSystem::phaseModelTable,  this->phaseModels_, iteri)
    {
        phaseModel& phasei = iteri()();

        phaseSystem::phaseModelTable::iterator iterk = iteri;
        iterk++;
        for
        (
            ;
            iterk != this->phaseModels_.end();
            ++iterk
        )
        {
            if (iteri()().name() != iterk()().name())
            {
                phaseModel& phasek = iterk()();

                // Phase i to phase k
                const phasePairKey keyik(phasei.name(), phasek.name(), true);

                // Phase k to phase i
                const phasePairKey keyki(phasek.name(), phasei.name(), true);

                // Net mass transfer from k to i phase
                tmp<volScalarField> tdmdtNetki
                (
                    new volScalarField
                    (
                        IOobject
                        (
                            "tdmdtYki",
                            this->mesh().time().timeName(),
                            this->mesh(),
                            IOobject::NO_READ,
                            IOobject::NO_WRITE
                        ),
                        this->mesh(),
                        dimensionedScalar
                        (
                            "zero",
                            dimDensity/dimTime,
                            0
                        )
                    )
                );
                volScalarField& dmdtNetki = tdmdtNetki.ref();


                if (massTransferModels_.found(keyik))
                {
                    autoPtr<interfaceCompositionModel>& interfacePtr =
                        massTransferModels_[keyik];

                    // Explicit temperature mass transfer rate
                    tmp<volScalarField> Kexp =
                        interfacePtr->Kexp
                        (
                            interfaceCompositionModel::T,
                            T
                        );

                    if (Kexp.valid())
                    {
                        Info << "Explicit temperature mass transfer.." << endl;
                        Info << "keyik :" << keyik << endl;

                        // Add explicit T based to all the other explixit terms
                        dmdtNetki -= Kexp.ref();
                        *dmdt_[keyik] = Kexp.ref();

                    }
                }

                // Looking for mass transfer in the other direction (k to i)
                if (massTransferModels_.found(keyki))
                {
                    autoPtr<interfaceCompositionModel>& interfacePtr =
                        massTransferModels_[keyki];

                    // Explicit temperature mass transfer rate
                    const tmp<volScalarField> Kexp =
                        interfacePtr->Kexp
                        (
                            interfaceCompositionModel::T,
                            T
                        );

                    if (Kexp.valid())
                    {
                        Info << "Explicit temperature mass transfer.." << endl;
                        Info << "keyki :" << keyki << endl;

                        dmdtNetki += Kexp.ref();
                        *dmdt_[keyki] = Kexp.ref();
                    }

                }

                word keyikName(phasei.name() + phasek.name());
                word keykiName(phasek.name() + phasei.name());

                eqn -=
                        (
                            dmdtNetki
                           *(
                                calculateL(dmdtNetki, keyik, keyki, T)
                              - (phasek.Cp() - phasei.Cp())
                              * dimensionedScalar("T0", dimTemperature, 298.0)
                            )
                        );
            }
        }
    }

    return tEqnPtr;
}


template<class BasePhaseSystem>
void Foam::MassTransferPhaseSystem<BasePhaseSystem>::massSpeciesTransfer
(
    const phaseModel& phase,
    volScalarField::Internal& Su,
    volScalarField::Internal& Sp,
    const word speciesName
)
{
    // Fill the volumetric mass transfer for species
    forAllIter(massTransferModelTable, massTransferModels_, iter)
    {
        if (iter()->transferSpecie() == speciesName)
        {
            // Extract alpha*div(u) added in alpha's
            tmp<volScalarField> divU = fvc::div(this->phi());

            // Explicit source
            Su =
                  this->Su()[phase.name()]
                - divU.ref().internalField()*phase.oldTime()
                + this->Sp()[phase.name()]*phase;

            // Implicit source
            //Sp = this->Sp()[phase.name()];
        }
    }
}

/*
template<class BasePhaseSystem>
void Foam::MassTransferPhaseSystem<BasePhaseSystem>::massTransfer
(
    const volScalarField& T
)
{
    forAllIter(massTransferModelTable, massTransferModels_, iter)
    {
        const phasePair& pair
        (
            this->phasePairs_[iter.key()]
        );

        // Disperse to to direction mass transfer
        const phaseModel& to = pair.to();
        const phaseModel& from = pair.from();

        const phasePairKey key(from.name(), to.name(), true);
        volScalarField& dmdt(*dmdt_[key]);

        // Explicit phase mass transfer rate
        tmp<volScalarField> Kexp =
            iter()->Kexp(interfaceCompositionModel::T, T);

        // Mass transfer with species
        if (iter()->transferSpecie() != "none")
        {
            if (Kexp.valid())
            {
                volScalarField& dmdtYiY =
                    dmdtYi
                    (
                        iter.key(),
                        interfaceCompositionModel::modelVariableNames
                        [
                            interfaceCompositionModel::T
                        ],
                        speciesName
                    );

                *dmdt_[key] = Kexp.ref();
                dmdtYiY = Kexp.ref();
            }
        }
        else
        {
            if (Kexp.valid())
            {
                *dmdt_[key] = Kexp.ref();
            }
        }
    }
}
*/
// ************************************************************************* //
