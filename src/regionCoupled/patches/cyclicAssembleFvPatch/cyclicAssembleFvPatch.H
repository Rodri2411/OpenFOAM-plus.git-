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
    Foam::cyclicAssembleFvPatch

Description
    Assemble fvPatch to cache nbrCells in global addressing for assemble matrix.
    Re-implementing internalFieldTransfer called from the GAMG ldu agglomeration
    function.

SourceFiles
    cyclicAssembleFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicAssembleFvPatch_H
#define cyclicAssembleFvPatch_H

#include "cyclicFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class cyclicAssembleFvPatch Declaration
\*---------------------------------------------------------------------------*/

class cyclicAssembleFvPatch
:
    public cyclicFvPatch
{
     // Private data

        //- Nbr faceCell
        const labelList& nrbFaceCells_;


public:

    //- Runtime type information
    TypeName(cyclicFvPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        cyclicAssembleFvPatch
        (
            const polyPatch& patch,
            const fvBoundaryMesh& bm,
            const labelList& nrbFaceCells
        )
        :
            cyclicFvPatch(patch, bm),
            nrbFaceCells_(nrbFaceCells)
        {}

    //- Destructor
    virtual ~cyclicAssembleFvPatch();


    // Member Functions


        //- Return neighbour field
        virtual tmp<labelField> internalFieldTransfer
        (
            const Pstream::commsTypes commsType,
            const labelUList& internalData
        ) const;



};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
