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
#include "volFields.H"

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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::assembleDecomp::findNbrMeshId
(
    const mappedPatchBase& pp,
    const UPtrList<fvMesh>& meshes
) const
{
    for(label i = 0; i < meshes.size(); i++)
    {
        if (meshes[i].name() == pp.sampleMesh().name())
        {
            return i;
        }
    }
    return -1;
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

    const wordList& fluidNames(rp["fluid"]);
    const wordList& solidsNames(rp["solid"]);

     //- Use regionProperties to find all other regions
    UPtrList<fvMesh> meshes(fluidNames.size() + solidsNames.size());

    wordList meshNames(meshes.size());

    label regioni = 0;
    forAll(fluidNames, i)
    {
        meshNames[regioni++] = fluidNames[i];
    }
    forAll(solidsNames, i)
    {
        meshNames[regioni++] = solidsNames[i];
    }

    forAll(meshes, i)
    {
        meshes.set
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

    lduPrimitiveMeshAssemble assembleLduMesh(meshes, io);

    const lduAddressing& addr = assembleLduMesh.lduAddr();

    globalIndex globalNumbering
    (
        addr.size(),
        Pstream::msgType(),
        mesh.comm(),
        Pstream::parRun()
    );

    const labelListList assembleCellCells
    (
        assembleLduMesh.globalCellCells
        (
            assembleLduMesh,
            globalNumbering
        )
    );

    vectorField cellCentres(addr.size(), vector::zero);

    forAll(meshes, i)
    {
        const label cellOffset = assembleLduMesh.cellOffsets()[i];

        forAll (meshes[i].C(), localCellI)
        {
            cellCentres[cellOffset + localCellI] = meshes[i].C()[localCellI];
        }
    }

    labelList assemblyDecomp
    (
        method_().decompose(assembleCellCells, cellCentres)
    );

    const labelUList& nbr = addr.upperAddr();
    const labelUList& own = addr.lowerAddr();

    PackedBoolList isMappedFace(addr.upperAddr().size());
    {
        forAll(meshes, i)
        {
            forAll(meshes[i].boundaryMesh(), patchI)
            {
                const polyPatch& pp = meshes[i].boundaryMesh()[patchI];
                if (isA<mappedPatchBase>(pp))
                {
                    forAll (pp, localFaceI)
                    {
                        label allFacei =
                            assembleLduMesh.faceBoundMap()[i][patchI]
                            [
                                localFaceI
                            ];
                        isMappedFace[allFacei] = true;
                    }
                }
            }
        }
    }

    while (true)
    {

        label nChanged = 0;

        forAll(own, facei)
        {
            if (isMappedFace[facei])
            {
                label ownProc = assemblyDecomp[own[facei]];
                label neiProc = assemblyDecomp[nbr[facei]];
                if (ownProc < neiProc)
                {
                    assemblyDecomp[nbr[facei]] = ownProc;
                    nChanged++;
                }
                else if (neiProc < ownProc)
                {
                    assemblyDecomp[own[facei]] = neiProc;
                    nChanged++;
                }
            }
        }

        if (returnReduce(nChanged, sumOp<label>()) == 0)
        {
            break;
        }
    }

    label offset(0);

    labelList localDecomp(mesh.nCells());

    forAll(meshes, i)
    {
        // Write decompositions for other regions
        if (mesh.name() != meshes[i].name())
        {
            labelIOList otherDecomp
            (
                IOobject
                (
                    dataFile_,
                    meshes[i].facesInstance(),
                    meshes[i],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                SubList<label>
                (
                    assemblyDecomp,
                    meshes[i].nCells(),
                    offset
                )
            );

            otherDecomp.write();
        }
        else
        {
            localDecomp = SubList<label>(assemblyDecomp, mesh.nCells(), offset);
        }

        offset += meshes[i].nCells();
    }

    // Return decomposition of my region
    return localDecomp;
}

// ************************************************************************* //
