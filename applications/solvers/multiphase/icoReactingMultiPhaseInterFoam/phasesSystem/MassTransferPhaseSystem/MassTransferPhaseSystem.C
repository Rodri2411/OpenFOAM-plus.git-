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

        L = -interfacePtr->L(species, T);

        return tL;
    }
    else if (massTransferModels_.found(keyki))
    {
        const autoPtr<interfaceCompositionModel>& interfacePtr =
            massTransferModels_[keyki];

        word speciesName = interfacePtr->transferSpecie();

        const label tempOpen(speciesName.find('.'));

        const word species(speciesName(0, tempOpen));

        L = interfacePtr->L(species, T);

        return tL;
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
//     tmp<fvScalarMatrix> tEqnPtr
//     (
//         new fvScalarMatrix(T, dimVolume*dimTemperature/dimTime)
//     );
    tmp<fvScalarMatrix> tEqnPtr
    (
        new fvScalarMatrix(T, dimEnergy/dimTime)
    );

    fvScalarMatrix& eqn = tEqnPtr.ref();

    forAllIter(phaseSystem::phaseModelTable,  this->phaseModels_, iteri)
    {
        phaseModel& phasei = iteri()();
//         const volScalarField invRhoCpi
//         (
//             1.0/phasei.rho()/phasei.Cp()
//         );

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

                //const volScalarField invRhoCpk(1.0/phasek.rho()/phasek.Cp());

                // Phase i to phase k
                const phasePairKey keyik(phasei.name(), phasek.name(), true);

                // Phase k to phase i
                const phasePairKey keyki(phasek.name(), phasei.name(), true);


                // mass transfer models coeffs*(T - Tref) for phase i
//                 tmp<volScalarField> tcoeffsi
//                 (
//                     new volScalarField
//                     (
//                         IOobject
//                         (
//                             "tcoeffsi",
//                             this->mesh().time().timeName(),
//                             this->mesh(),
//                             IOobject::NO_READ,
//                             IOobject::NO_WRITE
//                         ),
//                         this->mesh(),
//                         dimensionedScalar
//                         (
//                             "zero",
//                             dimDensity/dimTemperature/dimTime,
//                             0
//                         )
//                     )
//                 );
//                 volScalarField& coeffsi = tcoeffsi.ref();

                // mass transfer models coeffs*(T - Tref) for phase k
//                 tmp<volScalarField> tcoeffsk
//                 (
//                     new volScalarField
//                     (
//                         IOobject
//                         (
//                             "tcoeffsk",
//                             this->mesh().time().timeName(),
//                             this->mesh(),
//                             IOobject::NO_READ,
//                             IOobject::NO_WRITE
//                         ),
//                         this->mesh(),
//                         dimensionedScalar
//                         (
//                             "zero",
//                             dimDensity/dimTemperature/dimTime,
//                             0
//                         )
//                     )
//                 );
//                 volScalarField& coeffsk = tcoeffsk.ref();

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

                    //word speciesName = interfacePtr->transferSpecie();

                    //dSikdVar = interfacePtr->dSdVariable();

                    // Look for mass transfer i to k specie driven by Y and p
//                     if (dmdtYi_.found(keyik) || speciesName != "none")
//                     {
//                         Pair<word> YSpecie
//                         (
//                             interfaceCompositionModel::modelVariableNames
//                             [
//                                 interfaceCompositionModel::Y
//                             ],
//                             speciesName
//                         );
//
//                         if (dmdtYi_[keyik].found(YSpecie))
//                         {
//                             dmdtYki -= *dmdtYi_[keyik][YSpecie];
//                         }
//
//                         Pair<word> PSpecie
//                         (
//                             interfaceCompositionModel::modelVariableNames
//                             [
//                                 interfaceCompositionModel::P
//                             ],
//                             speciesName
//                         );
//
//                         if (dmdtYi_[keyik].found(PSpecie))
//                         {
//                             dmdtYki -= *dmdtYi_[keyik][PSpecie];
//                         }
//                     }


//                     volScalarField& dmdtYiT =
//                         dmdtYi
//                         (
//                             keyik,
//                             interfaceCompositionModel::modelVariableNames
//                             [
//                                 interfaceCompositionModel::T
//                             ],
//                             speciesName
//                         );

                    //dimensionedScalar Tref(interfacePtr->Tactivate());

                    // Fill dmdt for alpha's
                    //volScalarField& dmdtik = *dmdt_[keyik];

                    if (!interfacePtr->semiImplicit())
                    {
                        // Explicit temperature mass transfer rate
                        tmp<volScalarField> Kexp =
                            interfacePtr->Kexp
                            (
                                interfaceCompositionModel::T,
                                T
                            );

                        //if (Kexp.valid())
                        {
                            Info << "Explicit temperature mass transfer.." << endl;
                            Info << "keyik :" << keyik << endl;
                            // Add explicit T based to all the other explixit terms
                            dmdtNetki -= Kexp.ref();
                            //dmdtYki -= *dmdt_[keyik];
                            *dmdt_[keyik] = Kexp.ref();



                            // Add explicit source to dmdt_
                            //dmdtik = Kexp.ref();
                            // Add explicit source to dmdtYi to be used in YEq's
                            //dmdtYiT = Kexp.ref();

//                             Info<< " max(dmdtik) "
//                                 << max(dmdtik.internalField()) << endl;
//                             Info<< "temperature based Mass rate [kg/s]: "   << keyik
//                                 << gSum((dmdtik*this->mesh().V())()) << endl;
                        }
                    }
                    else if (interfacePtr->semiImplicit())
                    {
                        // Implicit temperature based mass transfer type :
                        // Kimp*(T - Tref)
                        // NOTE: dmdtYiT for species and dmdt for alpha's are
                        // lagged in T
                        tmp<volScalarField> Kimp =
                            interfacePtr->Kimp
                            (
                                interfaceCompositionModel::T,
                                T
                            );

                        if (Kimp.valid())
                        {
                            Info << "Implicit temperature mass transfer.." << endl;
                            Info << "keyik :" << keyik << endl;
                            // Fill alpha mass transfer
                            //dmdtik = Kimp.ref()*mag(T.oldTime() - Tref);
                            // Fill species mass transfer
                            //dmdtYiT = Kimp.ref()*mag(T.oldTime() - Tref);
                            //coeffsi = Kimp.ref();
                            ///Trefi = Tref;

//                             Info<< " max(dmdtik) "
//                                 << max(dmdtik.internalField()) << endl;
//                             Info<< "temperature based Mass rate [kg/s]: "   << keyik << " "
//                                 << gSum((dmdtik*this->mesh().V())()) << endl;
                        }
                    }
                }

                // derivate of source on variable (dSki/dVar) from k to i
                //label dSkidVar(0);

                // Looking for mass transfer in the other direction (k to i)
                if (massTransferModels_.found(keyki))
                {
                    autoPtr<interfaceCompositionModel>& interfacePtr =
                        massTransferModels_[keyki];

                    //word speciesName = interfacePtr->transferSpecie();

                    //dSkidVar = interfacePtr->dSdVariable();

                    // Look for mass transfer from k to i specie and pressure driven
//                     if (dmdtYi_.found(keyki) && speciesName != "none")
//                     {
//                         Pair<word> YSpecie
//                         (
//                             interfaceCompositionModel::modelVariableNames
//                             [
//                                 interfaceCompositionModel::Y
//                             ],
//                             speciesName
//                         );
//
//                         if (dmdtYi_[keyki].found(YSpecie))
//                         {
//                             dmdtYki += *dmdtYi_[keyki][YSpecie];
//                         }
//
//                         Pair<word> PSpecie
//                         (
//                             interfaceCompositionModel::modelVariableNames
//                             [
//                                 interfaceCompositionModel::P
//                             ],
//                             speciesName
//                         );
//
//                         if (dmdtYi_[keyki].found(PSpecie))
//                         {
//                             dmdtYki += *dmdtYi_[keyki][PSpecie];
//                         }
//                     }

                    //dimensionedScalar Tref(interfacePtr->Tactivate());
                    //volScalarField& dmdtki = *dmdt_(keyki);

//                     volScalarField& dmdtYiT =
//                         dmdtYi
//                         (
//                             keyki,
//                             interfaceCompositionModel::modelVariableNames
//                             [
//                                 interfaceCompositionModel::T
//                             ],
//                             speciesName
//                        );

                    if (!interfacePtr->semiImplicit())
                    {
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
                            //dmdtYki += *dmdt_[keyki];
                            *dmdt_[keyki] = Kexp.ref();


                            //dmdtki = Kexp.ref();
                            //dmdtYiT = Kexp.ref();

//                             Info<< " max(dmdtki) "
//                                 << max(dmdtki.internalField()) << endl;
                            Info<< " max(dmdtYki) "
                                << max(dmdtNetki.internalField()) << endl;
                            //Info<< " max(dmdtYiT) "
                            //    << max(dmdtYiT.internalField()) << endl;

//                             Info<< "temperature based Mass rate [kg/s]: "
//                                 << gSum((dmdtki*this->mesh().V())()) << endl;
                        }
                    }
                    else if (interfacePtr->semiImplicit())
                    {
                        // Implicit temperature based mass transfer
                        // NOTE: dmdtYiT and dmdt are lagged in T
                        tmp<volScalarField> Kimp =
                            interfacePtr->Kimp
                            (
                                interfaceCompositionModel::T,
                                T
                            );

                        if (Kimp.valid())
                        {
                            Info << "Implicit temperature mass transfer.." << endl;

                            //dmdtki = Kimp.ref()*mag(T.oldTime() - Tref);
                            //dmdtYiT = Kimp.ref()*mag(T.oldTime() - Tref);
                            //coeffsk = Kimp.ref();
                            //Trefk = Tref;

//                             Info<< " max(dmdtki) "
//                                 << max(dmdtki.internalField()) << endl;
//                             Info<< "temperature based Mass rate [kg/s]: "
//                                 << keyki << " "
//                                 << gSum((dmdtki*this->mesh().V())()) << endl;
                        }
                    }
                }
/*
                 eqn +=
                    (
                          dmdtYki
//                        + (-dSkidVar)*coeffsk*Trefk
//                        + (dSikdVar)*coeffsi*Trefi
                    )*tL.ref()
                    * (invRhoCpk - invRhoCpi)

                    - (
                        dmdtYki
//                      + (-dSkidVar)*coeffsk*Trefk
//                      + (dSikdVar)*coeffsi*Trefi
                      )*
                      (
                          invRhoCpk/iterk()->rho()
                        - invRhoCpi/iteri()->rho()
                      )*phasei.thermo().p();
*/
/*
                    + fvm::SuSp
                      (
                         ((dSkidVar)*coeffsk + (-dSikdVar)*coeffsi)*L
                        *(invRhoCpk - invRhoCpi)
                       - ((dSkidVar)*coeffsk + (-dSikdVar)*coeffsi)*
                         (
                            invRhoCpk/iterk()->rho()
                          - invRhoCpi/iteri()->rho()
                         )*phasei.thermo().p(),
                         T
                      );
*/

                eqn -=
                        (
                            dmdtNetki
                           *(
                              calculateL(keyik, keyki, T)
                            - (phasek.Cp() - phasei.Cp())
                              *dimensionedScalar("T0", dimTemperature, 298.0)
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
            //forAll (iter()->transferSpecie(), specieI)
            {
                //word speciesName = iter()->transferSpecie()[specieI];

                //const phasePair& pair(this->phasePairs_[iter.key()]);

                //const word from(pair.from().name());
                //const word to(pair.to().name());

                //const phasePairKey key(from, to, true);

                //const phaseModel& phaseFrom = this->phases()(from);
                //const volScalarField& phaseTo = phaseSystem::phases()(to);

//                 volScalarField& dmdtYiY =
//                     dmdtYi
//                     (
//                         key,
//                         interfaceCompositionModel::modelVariableNames
//                         [
//                             interfaceCompositionModel::T
//                         ],
//                         speciesName
//                     );
                DebugVar(max(this->Su()[phase.name()]));
                DebugVar(min(this->Su()[phase.name()]));

                DebugVar(max(this->Sp()[phase.name()]));
                DebugVar(min(this->Sp()[phase.name()]));

                // Extract alpha*div(u) added in alpha's
                tmp<volScalarField> divU = fvc::div(this->phi());

                // Explicit source
                Su =
                    this->Su()[phase.name()]
                  - divU.ref().internalField()*phase.oldTime()
                  - this->Sp()[phase.name()]*phase;

                // Implicit source
                //Sp = this->Sp()[phase.name()];//*phase/phase.oldTime();
            }
        }
    }
}


template<class BasePhaseSystem>
void Foam::MassTransferPhaseSystem<BasePhaseSystem>::massTransfer
(
    const volScalarField& T
)
{
    forAllConstIter(massTransferModelTable, massTransferModels_, iter)
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
            //if (!iter()->semiImplicit())
            {
                //forAll (iter()->transferSpecie(), specieI)
                {
                    //word speciesName = iter()->transferSpecie();

                    //word nameDisp(IOobject::groupName(speciesName, from.name()));
                    //word nameCont(IOobject::groupName(speciesName, to.name()));

                    //Info<< "key : " << key <<  " nameDisp to nameCont " << " "
                    //    << nameDisp << " " << nameCont << endl;

                    if (Kexp.valid())
                    {
//                         volScalarField& dmdtYiY =
//                             dmdtYi
//                             (
//                                 iter.key(),
//                                 interfaceCompositionModel::modelVariableNames
//                                 [
//                                     interfaceCompositionModel::T
//                                 ],
//                                 speciesName
//                             );

                        *dmdt_[key] = Kexp.ref();
                        //dmdtYiY = Kexp.ref();
                    }
                }
            }
        }
        else
        {
            if (Kexp.valid())
            {
                *dmdt_[key] = Kexp.ref();
            }
        }
            /*
            else if (iter()->semiImplicit())
            {
                // Implicit species based mass transfer
                tmp<volScalarField> Kimp = iter()->Kimp
                (
                    interfaceCompositionModel::Y,
                    T
                );

                if (Kimp.valid())
                {
                    Info << "Semi-implicit species based mass transfer.." << endl;

                    volScalarField& dmdtYiY =
                        dmdtYi
                        (
                            iter.key(),
                            interfaceCompositionModel::modelVariableNames
                            [
                                interfaceCompositionModel::Y
                            ],
                            speciesName
                        );

                    // Fill dm used in alpha's
                    dmdt = Kimp.ref()*
                    (
                        max
                        (
                            eqns[nameCont]->psi()
                          - iter()->Yf(speciesName, T),
                            0.0
                        )
                    );
                    dmdtYiY = dmdt;

                    *eqns[nameCont] +=
                      - Kimp.ref()*iter()->Yf(speciesName, T)
                      + fvm::Sp
                        (
                            Kimp.ref(),
                            eqns[nameCont]->psi()
                        );

                    // Explicit in the from phase
                    if (eqns.found(nameDisp))
                    {
                        *eqns[nameDisp] -= dmdt;
                    }
                }
            }
            */
        //}

        Info<< "Species based Mass rate [kg/s]: "
            << key << " "
            << gSum((dmdt*this->mesh().V())()) << endl;
    }
}

// ************************************************************************* //
