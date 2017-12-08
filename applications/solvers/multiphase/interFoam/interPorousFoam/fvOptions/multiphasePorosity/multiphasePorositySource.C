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

#include "multiphasePorositySource.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "fvmDdt.H"
#include "fvmSup.H"
#include "surfaceInterpolate.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(multiphasePorositySource, 0);
    addToRunTimeSelectionTable
    (
        option,
        multiphasePorositySource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::multiphasePorositySource::setPorosity()
{
    // Reset the porosiy field
    porosity_ == 1;

    forAll(zoneIDs_, i)
    {
        const labelList& zones = zoneIDs_[i];

        forAll(zones, j)
        {
            const label zonei = zones[j];
            const cellZone& cz = mesh_.cellZones()[zonei];

            forAll(cz, k)
            {
                const label celli = cz[k];
                porosity_[celli] = porosityZone_[i];
            }
        }
    }  

    porosity_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::multiphasePorositySource::multiphasePorositySource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(name, modelType, dict, mesh),
    aZone_(),
    bZone_(),
    cZone_(),
    D50Zone_(),
    porosityZone_(),
    zoneIDs_(),
    nuName_("nu"),
    porosityName_(coeffs_.lookup("porosity")),
    porosity_(mesh_.lookupObjectRef<volScalarField>(porosityName_)),
    kc_(-1)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::multiphasePorositySource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    const volVectorField& U = eqn.psi();
    const volScalarField& nu = mesh_.lookupObject<volScalarField>(nuName_);

    static bool writtenPorosity = false;
    if (!writtenPorosity && mesh_.time().writeTime())
    {
        porosity_.write();
        writtenPorosity = true;
    }

    volScalarField gamma
    (
        IOobject
        (
            typeName + ":gamma",
            mesh_.time().timeName(),
            mesh_.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("gamma", dimDensity/dimTime, 0)
    );

    volScalarField c
    (
        IOobject
        (
            typeName + ":c",
            mesh_.time().timeName(),
            mesh_.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("c", dimless, 0)
    );

    forAll(zoneIDs_, i)
    {
        const labelList& zones = zoneIDs_[i];

        forAll(zones, j)
        {
            const label zonei = zones[j];
            const cellZone& cz = mesh_.cellZones()[zonei];

            forAll(cz, k)
            {
                const label celli = cz[k];
                const scalar a = aZone_[i];
                const scalar b = bZone_[i];
                const scalar D50 = D50Zone_[i];

                const scalar p = porosity_[celli];
                const scalar rhoc = rho[celli];
                const vector& Uc = U[celli];
                const scalar nuc = nu[celli];

                scalar bTerm = b*mag(Uc);
                if (kc_ > 0)
                {
                    bTerm += 7.5/kc_*bTerm;
                }

                gamma[celli] =
                    (1.0 - p)*rhoc/D50/pow3(p)*(a*sqr(1.0 - p)*nuc/D50 + bTerm);

                c[celli] = cZone_[i];
            }
        }
    }  

    gamma.correctBoundaryConditions();
    c.correctBoundaryConditions();

    if (debug && mesh_.time().writeTime())
    {
        gamma.write();
        c.write();
    }

    fvMatrix<vector> porosityEqn
    (
      - c*fvm::ddt(rho, U)
      - fvm::Sp(porosity_*gamma, U)
    );

    // Contributions are added to RHS of momentum equation
    eqn += porosityEqn;
}


bool Foam::fv::multiphasePorositySource::read(const dictionary& dict)
{
    if (option::read(dict))
    {
        if (coeffs_.found("UNames"))
        {
            coeffs_.lookup("UNames") >> fieldNames_;
        }
        else if (coeffs_.found("U"))
        {
            word UName(coeffs_.lookup("U"));
            fieldNames_ = wordList(1, UName);
        }
        else
        {
            fieldNames_ = wordList(1, "U");
        }

        applied_.setSize(fieldNames_.size(), false);


        // Create the porosity models - 1 per region
        const dictionary& regionsDict(coeffs_.subDict("regions"));
        const wordList regionNames(regionsDict.toc());
        aZone_.setSize(regionNames.size(), 0);
        bZone_.setSize(regionNames.size(), 0);
        cZone_.setSize(regionNames.size(), 0);
        D50Zone_.setSize(regionNames.size(), 1);
        porosityZone_.setSize(regionNames.size(), 1);
        zoneIDs_.setSize(regionNames.size());
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
            modelDict.lookup("b") >> bZone_[i];
            modelDict.lookup("c") >> cZone_[i];
            modelDict.lookup("porosity") >> porosityZone_[i];
            modelDict.lookup("D50") >> D50Zone_[i];
        }

        // Field names
        // Note: porosity field name is invariant and only read on construction
        dict.readIfPresent("nu", nuName_);

        // Keulegan-Carpenter number
        if (dict.readIfPresent("kc", kc_))
        {
            Info<< "Employing Keulegan-Carpenter parameter, kc = " << kc_
                << endl;
        }

        // Set the porosity field
        setPorosity();

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
