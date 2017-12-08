/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
     \\/     M anipulation  | Copyright (C) 2017 IH-Cantabria
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

#include "mangrovesSource.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "fvmSup.H"
#include "fvmDdt.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(mangrovesSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        mangrovesSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::fv::mangrovesSource::kMangrove() const
{
    tmp<volScalarField> tk
    (
        new volScalarField
        (
            IOobject
            (
                name_ + ":k",
                mesh_.time().timeName(),
                mesh_.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("k", dimless/dimTime, 0)
        )
    );

    volScalarField& k = tk.ref();

    const volVectorField& U = mesh_.lookupObject<volVectorField>(UName_);

    forAll(zoneIDs_, i)
    {
        const labelList& zones = zoneIDs_[i];

        for (const label zonei : zones)
        {
            const cellZone& cz = mesh_.cellZones()[zonei];

            for (const label celli : cz)
            {
                const scalar a = aZone_[i];
                const scalar N = NZone_[i];
                const scalar Ckp = CkpZone_[i];
                const scalar Cd = CdZone_[i];

                const scalar magUc = mag(U[celli]);

                k[celli] = Ckp*Cd*a*N*magUc;
            }
        }
    }  

    k.correctBoundaryConditions();

    return tk;
}




Foam::tmp<Foam::volScalarField>
Foam::fv::mangrovesSource::epsilonMangrove() const
{
    tmp<volScalarField> tepsilon
    (
        new volScalarField
        (
            IOobject
            (
                name_ + ":epsilon",
                mesh_.time().timeName(),
                mesh_.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("epsilon", dimless/dimTime, 0)
        )
    );

    volScalarField& epsilon = tepsilon.ref();

    const volVectorField& U = mesh_.lookupObject<volVectorField>(UName_);

    forAll(zoneIDs_, i)
    {
        const labelList& zones = zoneIDs_[i];

        for (const label zonei : zones)
        {
            const cellZone& cz = mesh_.cellZones()[zonei];

            for (const label celli : cz)
            {
                const scalar a = aZone_[i];
                const scalar N = NZone_[i];
                const scalar Cep = CepZone_[i];
                const scalar Cd = CdZone_[i];

                const scalar magUc = mag(U[celli]);

                epsilon[celli] = Cep*Cd*a*N*magUc;
            }
        }
    }  

    epsilon.correctBoundaryConditions();

    return tepsilon;
}


Foam::tmp<Foam::volScalarField>
Foam::fv::mangrovesSource::UDragMangrove() const
{
    tmp<volScalarField> tUDrag
    (
        new volScalarField
        (
            IOobject
            (
                name_ + ":UDrag",
                mesh_.time().timeName(),
                mesh_.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("UDrag", dimless/dimTime, 0)
        )
    );

    volScalarField& UDrag = tUDrag.ref();

    const volVectorField& U = mesh_.lookupObject<volVectorField>(UName_);

    forAll(zoneIDs_, i)
    {
        const labelList& zones = zoneIDs_[i];

        for (const label zonei : zones)
        {
            const cellZone& cz = mesh_.cellZones()[zonei];

            for (const label celli : cz)
            {
                const scalar a = aZone_[i];
                const scalar N = NZone_[i];
                const scalar Cd = CdZone_[i];

                const scalar magUc = mag(U[celli]);

                UDrag[celli] = 0.5*Cd*a*N*magUc;
            }
        }
    }  

    UDrag.correctBoundaryConditions();

    return tUDrag;
}


Foam::tmp<Foam::volScalarField>
Foam::fv::mangrovesSource::UInertiaMangrove() const
{
    tmp<volScalarField> tUInertia
    (
        new volScalarField
        (
            IOobject
            (
                name_ + ":UInertia",
                mesh_.time().timeName(),
                mesh_.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("UInertia", dimless, 0)
        )
    );

    volScalarField& UInertia = tUInertia.ref();

    const scalar pi = constant::mathematical::pi;

    forAll(zoneIDs_, i)
    {
        const labelList& zones = zoneIDs_[i];

        for (const label zonei : zones)
        {
            const cellZone& cz = mesh_.cellZones()[zonei];

            for (const label celli : cz)
            {
                const scalar a = aZone_[i];
                const scalar Cm = CmZone_[i];

                UInertia[celli] = 0.25*(Cm + 1)*pi*sqr(a);
            }
        }
    }  

    UInertia.correctBoundaryConditions();

    return tUInertia;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::mangrovesSource::mangrovesSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(name, modelType, dict, mesh),
    aZone_(),
    NZone_(),
    CkpZone_(),
    CepZone_(),
    CmZone_(),
    CdZone_(),
    UName_("U"),
    kName_("k"),
    epsilonName_("epsilon")
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::mangrovesSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (eqn.psi().name() == epsilonName_)
    {
        DebugInFunction << "Applying epsilon source to field "
            << epsilonName_ << nl << endl;

        eqn -= fvm::Sp(epsilonMangrove(), eqn.psi());
    }
    else if (eqn.psi().name() == kName_)
    {
        DebugInFunction << "Applying k source to field "
            << kName_ << nl << endl;

        eqn -= fvm::Sp(kMangrove(), eqn.psi());
    }
}


void Foam::fv::mangrovesSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (eqn.psi().name() == epsilonName_)
    {
        DebugInFunction << "Applying epsilon source to field "
            << epsilonName_ << nl << endl;

        eqn -= fvm::Sp(rho*epsilonMangrove(), eqn.psi());
    }
    else if (eqn.psi().name() == kName_)
    {
        DebugInFunction
            << "Applying k source to field "
            << kName_ << nl << endl;

        eqn -= fvm::Sp(rho*kMangrove(), eqn.psi());
    }
}


void Foam::fv::mangrovesSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (eqn.psi().name() == epsilonName_)
    {
        DebugInFunction
            << "Applying epsilon source to field "
            << epsilonName_ << nl << endl;

        eqn -= fvm::Sp(alpha*rho*epsilonMangrove(), eqn.psi());
    }
    else if (eqn.psi().name() == kName_)
    {
        DebugInFunction
            << "Applying k source to field "
            << kName_ << nl << endl;

        eqn -= fvm::Sp(alpha*rho*kMangrove(), eqn.psi());
    }
}


void Foam::fv::mangrovesSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    DebugInFunction
        << "Applying momentum source to field "
        << eqn.psi().name() << nl << endl;

    eqn -=
        fvm::Sp(UDragMangrove(), eqn.psi())
      + UInertiaMangrove()*fvm::ddt(eqn.psi());
}


void Foam::fv::mangrovesSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    DebugInFunction
        << "Applying momentum source to field "
        << eqn.psi().name() << nl << endl;

    eqn -=
        fvm::Sp(rho*UDragMangrove(), eqn.psi())
      + UInertiaMangrove()*fvm::ddt(rho, eqn.psi());
}


void Foam::fv::mangrovesSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    DebugInFunction
        << "Applying momentum source to field "
        << eqn.psi().name() << nl << endl;

    eqn -=
        fvm::Sp(alpha*rho*UDragMangrove(), eqn.psi())
      + UInertiaMangrove()*fvm::ddt(alpha*rho, eqn.psi());
}


bool Foam::fv::mangrovesSource::read(const dictionary& dict)
{
    if (option::read(dict))
    {
        const dictionary& regionsDict(coeffs_.subDict("regions"));
        const wordList regionNames(regionsDict.toc());
        const label nRegion = regionNames.size();
        zoneIDs_.setSize(nRegion);

        bool turbulence = readBool(coeffs_.lookup("turbulence"));

        coeffs_.readIfPresent("U", UName_);
        if (turbulence)
        {
            coeffs_.readIfPresent("k", kName_);
            coeffs_.readIfPresent("epsilon", epsilonName_);

            fieldNames_.setSize(3);
            fieldNames_[0] = UName_;
            fieldNames_[1] = epsilonName_;
            fieldNames_[2] = kName_;

            CkpZone_.setSize(nRegion, 1);
            CepZone_.setSize(nRegion, 1);
        }
        else
        {
            fieldNames_.setSize(1);
            fieldNames_[0] = UName_;
        }

        applied_.setSize(fieldNames_.size(), false);

        // Create the Mangroves models - 1 per region
        aZone_.setSize(nRegion, 1);
        NZone_.setSize(nRegion, 1);
        CmZone_.setSize(nRegion, 1);
        CdZone_.setSize(nRegion, 1);

        forAll(zoneIDs_, i)
        {
            const word& regionName = regionNames[i];
            const dictionary& modelDict = regionsDict.subDict(regionName);

            const word& zoneName = modelDict.lookup("cellZone");
            zoneIDs_[i] = mesh_.cellZones().findIndices(zoneName);

            if (zoneIDs_[i].empty())
            {
                FatalErrorInFunction
                    << "Unable to find cellZone " << zoneName << nl
                    << "Valid cellZones are:" << mesh_.cellZones().names()
                    << exit(FatalError);
            }

            modelDict.lookup("a") >> aZone_[i];
            modelDict.lookup("N") >> NZone_[i];
            modelDict.lookup("Cm") >> CmZone_[i];
            modelDict.lookup("Cd") >> CdZone_[i];

            if (turbulence)
            {
                modelDict.lookup("Ckp") >> CkpZone_[i];
                modelDict.lookup("Cep") >> CepZone_[i];
            }
        }

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
