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

Class
    Foam::lduPrimitiveMeshAssembleAssemble

Description
    Assemble lduMatrix specific use for inter-region coupling through mapped
    patches.

SourceFiles
    lduPrimitiveMeshAssemble.C

\*---------------------------------------------------------------------------*/

#ifndef lduPrimitiveMeshAssemble_H
#define lduPrimitiveMeshAssemble_H

#include "fvMesh.H"
#include "lduPrimitiveMesh.H"
#include "mappedPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                       Class lduPrimitiveMeshAssemble Declaration
\*---------------------------------------------------------------------------*/

class lduPrimitiveMeshAssemble
:
    public objectRegistry,
    public lduPrimitiveMesh
{
    // Private data


        List<labelList> patchAddr_;

        //- Map patches from local to global
        labelListList patchMap_;

        //- Map faces from local to global
        labelListList faceMap_;

        //- Map path face from local to global
        labelListListList faceBoundMap_;

        // Off set cells
        labelList cellOffsets_;


    // Private Member Functions

        //- Get size of all meshes
        static label totalSize(const UPtrList<fvMesh>&);

        //- Find nrb mesh Id for mapped patches
        label findNbrMeshId
        (
            const mappedPatchBase&,
            const UPtrList<fvMesh>&
        ) const;

        //- Disallow default bitwise copy construct
        lduPrimitiveMeshAssemble
        (
            const lduPrimitiveMeshAssemble&,
            const IOobject& io
        );

        //- Disallow default bitwise assignment
        void operator=(const lduPrimitiveMeshAssemble&);


public:

    // Static data

        // Declare name of the class and its debug switch
        ClassName("lduPrimitiveMeshAssemble");

    // Constructors


        //- Construct from fvMeshes combining mapped faces into internal faces
        lduPrimitiveMeshAssemble
        (
            const UPtrList<fvMesh>& meshes,  const IOobject& io
        );


    //- Destructor
    virtual ~lduPrimitiveMeshAssemble()
    {}


    // Member Functions

        // Access

            //- Return patch addressing
            virtual const labelUList& patchAddr(const label patchNo) const
            {
                return patchAddr_[patchNo];
            }

             //- Return the object registry
            virtual const objectRegistry& thisDb() const
            {
                return *this;
            }

            //- Return patchMap
            const labelListList& patchMap() const
            {
                return patchMap_;
            }

            //- Return boundary face map
            const labelListListList& faceBoundMap() const
            {
                return faceBoundMap_;
            }

            //- Return cellOffsets
            const labelList& cellOffsets() const
            {
                return cellOffsets_;
            }

            //- Return faceMap
            const labelListList& faceMap() const
            {
                return faceMap_;
            }

            //- Return patchAddr
            const List<labelList>& patchAddr() const
            {
                return patchAddr_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
