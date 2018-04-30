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

#include "lduPrimitiveMeshAssemble.H"
#include "fvMesh.H"
#include "cyclicLduInterface.H"
#include "cyclicAssembleFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(lduPrimitiveMeshAssemble, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


Foam::label Foam::lduPrimitiveMeshAssemble::totalSize
(
    const UPtrList<fvMesh>& meshes
)
{
    label size = 0;

    forAll(meshes, i)
    {
        size += meshes[i].lduAddr().size();
    }
    return size;
}


Foam::label Foam::lduPrimitiveMeshAssemble::findNbrMeshId
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


Foam::lduPrimitiveMeshAssemble::lduPrimitiveMeshAssemble
(
    const UPtrList<fvMesh>& meshes,
    const IOobject& io
)
:
    objectRegistry(io),
    lduPrimitiveMesh(totalSize(meshes))
{
    forAll(meshes, i)
    {
        if (meshes[i].comm() != comm())
        {
            WarningInFunction
                << "Communicator " << meshes[i].comm()
                << " at index " << i
                << " differs between meshes "
                << endl;
        }
    }

    const label nMeshes = meshes.size();
    patchMap_.setSize(nMeshes);
    faceBoundMap_.setSize(nMeshes);
    faceMap_.setSize(nMeshes);

    // Determine cellOffset and faceOffset
    cellOffsets_.setSize(nMeshes);
    cellOffsets_[0] = 0;
    for(label i=1; i < nMeshes; i++)
    {
        cellOffsets_[i] = cellOffsets_[i-1] + meshes[i-1].lduAddr().size();
    }

    // Get newFaces and newPatches (not mapPolyPatches)
    label newFaces(0);
    label newPatches(0);
    for(label i=0; i < nMeshes; i++)
    {
        patchMap_[i].setSize(meshes[i].boundaryMesh().size(), -1);

        forAll (meshes[i].boundaryMesh(), patchI)
        {
            const polyPatch& pp = meshes[i].boundaryMesh()[patchI];
            if (!isA<mappedPatchBase>(pp))
            {
                patchMap_[i][patchI] = newPatches++;
            }
            else
            {
                const mappedPatchBase& mpp = refCast<const mappedPatchBase>(pp);

                const label meshNrbID = findNbrMeshId(mpp, meshes);

                const label nbrPatchID =
                    meshes[meshNrbID].boundaryMesh().findPatchID
                    (
                        mpp.samplePatch()
                    );

                const polyPatch& nbrpp =
                    meshes[meshNrbID].boundaryMesh()[nbrPatchID];

                if (pp.size() != nbrpp.size())
                {
                    FatalErrorInFunction
                        << "The number of faces on either side of the mapped"
                        << "patch " << pp.name() << " are not the same."
                        << "This might be due to the decomposition used.Please"
                        << " use the type assembleDecomp."
                        << exit(FatalError);
                }
                if (mpp.owner())
                {
                    newFaces += pp.size();
                }
            }
        }
    }

    // Add the internal faces for each mesh
    for(label i=0; i < nMeshes; i++)
    {
        newFaces += meshes[i].lduAddr().upperAddr().size();
    }

     // This gives the global cellId given the local patchId for interfaces
    patchAddr_.setSize(newPatches);

    for(label i=0; i < nMeshes; i++)
    {
        const lduInterfacePtrsList interfacesLst = meshes[i].interfaces();
        forAll(interfacesLst, patchI)
        {
            label globalPatchId = patchMap_[i][patchI];
            if (globalPatchId != -1)
            {
                const labelUList& faceCells =
                    meshes[i].lduAddr().patchAddr(patchI);

                // Fill local patchAddr for standard patches
                if (faceCells.size() > 0)
                {
                    patchAddr_[globalPatchId].setSize(faceCells.size(), -1);

                    for (label celli = 0; celli < faceCells.size(); celli++)
                    {
                        patchAddr_[globalPatchId][celli] =
                            cellOffsets_[i] + faceCells[celli];
                    }
                }
            }
        }
    }

    interfaces().setSize(newPatches);
    // Primitive interfaces
    primitiveInterfaces().setSize(newPatches);


    // The interfaces are conserved (cyclics, proc, etc)
    label interfaceID = 0;
    for(label i=0; i < nMeshes; i++)
    {
        const lduInterfacePtrsList interfacesLst = meshes[i].interfaces();

        faceBoundMap_[i].setSize(interfacesLst.size());

        forAll(interfacesLst, patchI)
        {
            label globalPatchId = patchMap_[i][patchI];
            if (globalPatchId != -1)
            {
                // Set interfaces for cyclic (cyclicAssembleFvPatch).
                // cyclic patches are cloned resetting nrbPatchID to the new
                // address plus local and nbr faceCells in global addressing
                // NOTE: the fvBoundary used remains the local. This means
                // that some member functions of the new cloned patch are
                // not fully functional, i.e neighbPatch() uses boundaryMesh
                // to access nbr patch.
                if (interfacesLst.set(patchI))
                {
                    if (isA<cyclicLduInterface>(interfacesLst[patchI]))
                    {
                        label nbrId = refCast
                            <const cyclicLduInterface>
                            (
                                interfacesLst[patchI]
                            ).neighbPatchID();

                        label globalNbr = patchMap()[i][nbrId];

                        primitiveInterfaces().set
                        (
                            interfaceID,
                            new cyclicAssembleFvPatch
                            (
                               *(
                                    new cyclicPolyPatch
                                    (
                                        refCast<const cyclicPolyPatch>
                                        (
                                            meshes[i].boundaryMesh()[patchI]
                                        ),
                                        globalNbr,
                                        patchAddr_[globalPatchId]
                                    )
                                ),
                                meshes[i].boundary(),
                                patchAddr_[globalNbr]
                            )
                        );

                        interfaces().set
                        (
                            interfaceID,
                            &primitiveInterfaces()[interfaceID]
                        );
                    }
                    else
                    {
                        primitiveInterfaces().set
                        (
                            interfaceID,
                            nullptr
                        );

                        interfaces().set
                        (
                            interfaceID,
                            interfacesLst(patchI)
                        );
                    }
                }
                interfaceID++;
            }
        }
    }

    // Create new addressing
    lowerAddr().setSize(newFaces, -1);
    upperAddr().setSize(newFaces, -1);

    label startIndex = 0;

    for(label i=0; i < nMeshes; i++)
    {
        faceMap_[i].setSize(meshes[i].lduAddr().lowerAddr().size(), -1);

        label localFacei = 0;

        label nFaces = meshes[i].lduAddr().upperAddr().size();

        // Add indivual addresses
        SubList<label>(lowerAddr(), nFaces, startIndex) =
            meshes[i].lduAddr().lowerAddr();

        SubList<label>(upperAddr(), nFaces, startIndex) =
            meshes[i].lduAddr().upperAddr();

        // Offset cellsID's to global cell addressing
        for (label facei=startIndex; facei < startIndex + nFaces; facei++)
        {
            lowerAddr()[facei] += cellOffsets_[i];
            upperAddr()[facei] += cellOffsets_[i];

            faceMap_[i][localFacei++] = facei;
        }

        startIndex += nFaces;
    }

    // Add lower/upper adressing for new internal faces corresponding
    // to old mapPolyPatch's
    label nFaces = startIndex;

    for(label i=0; i < nMeshes; i++)
    {
        forAll (meshes[i].boundaryMesh(), patchI)
        {
            const polyPatch& pp = meshes[i].boundaryMesh()[patchI];
            if (isA<mappedPatchBase>(pp))
            {
                const mappedPatchBase& mpp = refCast<const mappedPatchBase>(pp);

                if (mpp.owner())
                {
                    labelList nbrFaceCells(mpp.samplePolyPatch().faceCells());

                    const label meshNrbId = findNbrMeshId(mpp, meshes);

                    if (meshNrbId != -1)
                    {
                        forAll (pp.faceCells(), faceI)
                        {
                            label cellI = pp.faceCells()[faceI] + cellOffsets_[i];

                            label nbrCellI =
                                nbrFaceCells[faceI] + cellOffsets_[meshNrbId];

                            lowerAddr()[nFaces] = min(cellI, nbrCellI);
                            upperAddr()[nFaces] = max(cellI, nbrCellI);

                            nFaces++;
                        }
                    }
                    else
                    {
                         FatalErrorInFunction
                            << "Can not find Nbr Mesh for mapped patch"
                            << exit(FatalError);
                    }
                }
            }
        }
    }

    if (newFaces != nFaces)
    {
       FatalErrorInFunction
            << "The number of total faces in the assembled ldu matrix is wrong "
            << exit(FatalError);
    }

    // Fill faceBoundMap
    nFaces = startIndex;
    for(label i=0; i < nMeshes; i++)
    {
        forAll (meshes[i].boundaryMesh(), patchI)
        {
            const polyPatch& pp = meshes[i].boundaryMesh()[patchI];
            if (isA<mappedPatchBase>(pp) && pp.size() > 0)
            {
                const mappedPatchBase& mpp = refCast<const mappedPatchBase>(pp);
                if (mpp.owner())
                {
                    const label samplePatchi = mpp.samplePolyPatch().index();
                    label meshNrbId = findNbrMeshId(mpp, meshes);
                    const polyPatch& ppNbr = mpp.samplePolyPatch();

                    if
                    (
                        faceBoundMap_[i][patchI].size() == 0
                     && faceBoundMap_[meshNrbId][samplePatchi].size() == 0
                    )
                    {
                        faceBoundMap_[i][patchI].setSize(pp.size(), -1);
                        faceBoundMap_[meshNrbId][samplePatchi].setSize
                        (
                            ppNbr.size(),
                            -1
                        );

                        forAll (pp.faceCells(), faceI)
                        {
                            faceBoundMap_[i][patchI][faceI] = nFaces;
                            faceBoundMap_[meshNrbId][samplePatchi][faceI] = nFaces;
                            nFaces++;
                        }
                    }
                }
            }
        }
    }

    // Sort upper-tri order
    {
        labelList oldToNew
        (
            upperTriOrder
            (
                lduAddr().size(),
                lowerAddr(),
                upperAddr()
            )
        );

        inplaceReorder(oldToNew, lowerAddr());
        inplaceReorder(oldToNew, upperAddr());

        forAll(faceBoundMap_, meshi)
        {
            labelListList& bMap = faceBoundMap_[meshi];
            forAll(bMap, bi)
            {
                labelList& faceMap = bMap[bi];
                forAll(faceMap, i)
                {
                    faceMap[i] = oldToNew[faceMap[i]];
                }
            }
        }

        forAll(faceMap_, meshi)
        {
            labelList& faceMap = faceMap_[meshi];
            forAll(faceMap, facei)
            {
                faceMap[facei] = oldToNew[faceMap[facei]];
            }
        }
    }

    if (debug)
    {
        //DebugVar(faceBoundMap_);
        //DebugVar(lowerAddr());
        //DebugVar(upperAddr());
        //DebugVar(patchAddr_);
        //DebugVar(cellOffsets_);
        //DebugVar(faceMap_);
        checkUpperTriangular(lduAddr().size(), lowerAddr(), upperAddr());
    }
}

// ************************************************************************* //

