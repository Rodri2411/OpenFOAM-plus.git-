/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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

#include "assembleDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "globalIndex.H"
#include "regionProperties.H"
#include "fvMesh.H"
#include "lduPrimitiveMeshAssemble.H"
#include "labelIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(assembleDecomp, 0);

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        assembleDecomp,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::assembleDecomp::assembleDecomp(const dictionary& decompDict)
:
    decompositionMethod(decompDict),
    methodDict_(findCoeffsDict(typeName + "Coeffs", selectionType::MANDATORY)),
    dataFile_(findCoeffsDict(typeName + "Coeffs").lookup("dataFile"))
{
    method_ = decompositionMethod::New(methodDict_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::assembleDecomp::decompose
(
    const polyMesh& mesh,
    const pointField& points,
    const scalarField& pointWeights
)
{

    regionProperties rp(mesh.time());

    const wordList fluidNames(rp["fluid"]);
    const wordList solidsNames(rp["solid"]);

     //- Use regionProperties to find all other regions
    UPtrList<fvMesh> allMeshes(fluidNames.size() + solidsNames.size());

    wordList meshNames(allMeshes.size());

    label regioni = 0;
    forAll(fluidNames, i)
    {
        meshNames[regioni++] = fluidNames[i];
    }
    forAll(solidsNames, i)
    {
        meshNames[regioni++] = solidsNames[i];
    }

    forAll(allMeshes, i)
    {
        allMeshes.set
        (
            i,
            new fvMesh
            (
                IOobject
                (
                    meshNames[i],
                    mesh.time().timeName(),
                    mesh.time(),
                    IOobject::MUST_READ
                )
            )
        );
    }

    IOobject io
    (
        "assembleLdu",
        mesh.time().timeName(),
        mesh.time(),
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
    );

    lduPrimitiveMeshAssemble assembleLduMesh(allMeshes, io);

    const lduAddressing& addr = assembleLduMesh.lduAddr();

    globalIndex globalNumbering
    (
        addr.size(),
        Pstream::msgType(),
        mesh.comm(),
        Pstream::parRun()
    );

    labelListList assembleCellCells =
        static_cast<const lduPrimitiveMesh&>(assembleLduMesh).globalCellCells
        (
            assembleLduMesh,
            globalNumbering
        );

    const pointField dummyCc(addr.size(), vector::zero);

    labelList subDecomp
    (
        method_().decompose(assembleCellCells, dummyCc)
    );

    const labelUList& nbr = addr.upperAddr();
    const labelUList& own = addr.lowerAddr();

    forAll(allMeshes, i)
    {
        // extend the same processor Id on each side of the mappedPatch
        forAll(allMeshes[i].boundaryMesh(), patchI)
        {
            const polyPatch& pp = allMeshes[i].boundaryMesh()[patchI];
            if (isA<mappedPatchBase>(pp))
            {
                const mappedPatchBase& mpp =
                    refCast<const mappedPatchBase>(pp);

                if (mpp.owner())
                {
                    forAll (pp, localFaceI)
                    {
                        //label faceI = pp.start() + localFaceI;
                        label allFacei =
                            assembleLduMesh.faceBoundMap()[i][patchI]
                            [
                                localFaceI
                            ];

                        subDecomp[nbr[allFacei]] = subDecomp[own[allFacei]];
                    }
                }
            }
        }
    }

    label offset(0);

    forAll(allMeshes, i)
    {
        // Write decompositions for other regions
        if (mesh.name() != allMeshes[i].name())
        {
            labelIOList otherDecomp
            (
                IOobject
                (
                    dataFile_,
                    allMeshes[i].facesInstance(),
                    allMeshes[i],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                SubList<label>
                (
                    subDecomp,
                    allMeshes[i].nCells(),
                    offset
                )
            );

            otherDecomp.write();
        }

        offset += allMeshes[i].nCells();
    }

    const labelList localDecomp = SubList<label>(subDecomp, mesh.nCells());

    // Return decomposition of my region
    return localDecomp;
}

// ************************************************************************* //
