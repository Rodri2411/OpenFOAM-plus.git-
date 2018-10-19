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

#include "surfaceWriterBase.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::pointField Foam::surfaceWriters::writerBase::dummyPoints_;
Foam::faceList Foam::surfaceWriters::writerBase::dummyFaces_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceWriters::writerBase::writerBase()
:
    writerBase(dummyPoints_, dummyFaces_)
{}


Foam::surfaceWriters::writerBase::writerBase
(
    const dictionary& options
)
:
    writerBase()
{}


Foam::surfaceWriters::writerBase::writerBase
(
    const pointField& points,
    const faceList& faces
)
:
    pointsRef_(std::cref<pointField>(points)),
    facesRef_(std::cref<faceList>(faces))
{}


Foam::surfaceWriters::writerBase::writerBase
(
    const pointField& points,
    const faceList& faces,
    const dictionary& options
)
:
    writerBase(points, faces)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// void Foam::vtk::surfaceWriter::setTime(const instant& inst)
// {
//     instant_ = inst;
// }
//
//
// bool Foam::vtk::surfaceWriter::beginFile(std::string title)
// {
//     if (title.size())
//     {
//         return vtk::fileWriter::beginFile(title);
//     }
//
//     if (instant_.name().size())
//     {
//         return vtk::fileWriter::beginFile
//         (
//             "time='" + instant_.name() + "'"
//         );
//     }
//
//     // Provide default title
//     return vtk::fileWriter::beginFile("surface");
// }
//
//
// bool Foam::vtk::surfaceWriter::writeGeometry()
// {
//     enter_Piece();
//
//     beginPiece();
//
//     writePoints();
//
//     const globalIndex globalPointOffset(nLocalPoints_);
//
//     if (legacy())
//     {
//         writePolysLegacy(globalPointOffset);
//     }
//     else
//     {
//         writePolys(globalPointOffset);
//     }
//
//     return true;
// }
//
//
// bool Foam::vtk::surfaceWriter::beginCellData(label nFields)
// {
//     return enter_CellData(numberOfCells_, nFields);
// }
//
//
// bool Foam::vtk::surfaceWriter::beginPointData(label nFields)
// {
//     return enter_PointData(numberOfPoints_, nFields);
// }
//
//
// void Foam::vtk::surfaceWriter::writeTimeValue()
// {
//     if (instant_.name().size())
//     {
//         vtk::fileWriter::writeTimeValue(instant_.value());
//     }
// }


// ************************************************************************* //
