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
    Foam::vtk::surfaceWriter

Description
    Write faces/points (optionally with fields)
    as a vtp file or a legacy vtk file.

    The file output states are managed by the Foam::vtk::fileWriter class.
    FieldData (eg, TimeValue) must appear before any geometry pieces.

Note
    Parallel output is combined into a single Piece without point merging,
    which is similar to using multi-piece data sets, but allows more
    convenient creation as a streaming process.
    In the future, the duplicate points at processor connections
    may be addressed using ghost points.

SourceFiles
    foamVtkSurfaceWriter.C
    foamVtkSurfaceWriterTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceWriterBase_H
#define surfaceWriterBase_H

#include <functional>
#include "foamVtkFileWriter.H"
#include "pointField.H"
#include "faceList.H"
#include "instant.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class globalIndex;

namespace surfaceWriters
{

/*---------------------------------------------------------------------------*\
                     Class vtk::surfaceWriter Declaration
\*---------------------------------------------------------------------------*/

class writerBase
{
protected:

    // Protected data

        static pointField dummyPoints_;
        static faceList dummyFaces_;

        //- Reference to the points
        std::reference_wrapper<const pointField> pointsRef_;

        //- Reference to the faces
        std::reference_wrapper<const faceList> facesRef_;


    // Private Member Functions

        //- Determing sizes (nLocalPoints_, nLocalFaces_, nLocalVerts_),
        //- and begin piece
        void beginPiece();

        //- Write patch points
        void writePoints();

        //- Write patch faces, legacy format
        void writePolysLegacy(const globalIndex& pointOffsets);

        //- Write patch faces
        void writePolys(const globalIndex& pointOffsets);


public:

    // Constructors

        //- Construct from components (default format INLINE_BASE64)
        writerBase();

        //- Construct from components and open the file for writing.
        //  The file name is with/without an extension.
        writerBase(const dictionary& options);

        //- Construct from components (default format INLINE_BASE64),
        //- and open the file for writing.
        //  The file name is with/without an extension.
        writerBase
        (
            const pointField& points,
            const faceList& faces
        );

        //- Construct from components and open the file for writing.
        //  The file name is with/without an extension.
        writerBase
        (
            const pointField& points,
            const faceList& faces,
            const dictionary& options
        );


    //- Destructor
    virtual ~writerBase() = default;


    // Member Functions

        //- Define a time name/value for the output
        virtual void setTime(const instant& inst)     { };

        //- Write file header (non-collective)
        //  \note Expected calling states: (OPENED).
        virtual bool beginFile(std::string title = "")   { return true; };

        //- Write patch topology
        //  Also writes the file header if not previously written.
        //  \note Must be called prior to writing CellData or PointData
        virtual bool writeGeometry()   { return true; };

        //- Begin CellData output section for specified number of fields.
        //  Must be called prior to writing any cell data fields.
        //  \param nFields is for legacy format only.
        //      When nFields=0, this a no-op for legacy format.
        //  \note Expected calling states: (PIECE | POINT_DATA).
        //
        //  \return True if the state changed
        virtual bool beginCellData(label nFields = 0)   { return true; };

        //- Begin PointData for specified number of fields.
        //  Must be called prior to writing any point data fields.
        //  \param nFields is for legacy format only.
        //      When nFields=0, this a no-op for legacy format.
        //  \note Expected calling states: (PIECE | CELL_DATA).
        //
        //  \return True if the state changed
        virtual bool beginPointData(label nFields = 0)   { return true; };


        // Write
//
//         //- Write a list of Cell (Face) or Point values
//         template<class Type>
//         void write(const word& fieldName, const UList<Type>& field);

//         void write(const word& fieldName, const UList<scalar>& field);
//         void write(const word& fieldName, const UList<vector>& field);
//         void write(const word& fieldName, const UList<sphericalTensor>& field);
//         void write(const word& fieldName, const UList<symmTensor>& field);
//         void write(const word& fieldName, const UList<tensor>& field);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceWriters
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
