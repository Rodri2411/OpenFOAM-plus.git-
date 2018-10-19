/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2018 OpenCFD Ltd.
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

#include "rawSurfaceWriterSeq.H"
#include "makeSurfaceWriterMethods.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceWriters
{
     defineTypeNameAndDebug(rawSeq, 0);
//     defineRunTimeSelectionTable(factory, word);
//     defineRunTimeSelectionTable(factory, wordDict);
//     addNamedToRunTimeSelectionTable
//     (
//         factory,
//         factory,
//         word,
//         null
//     );

}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

#if 0
namespace Foam
{
    template<>
    void Foam::rawSurfaceWriter::writeHeader
    (
        Ostream& os,
        const word& fieldName,
        const Field<scalar>& values
    )
    {
        os  << values.size() << nl
            << "#  x  y  z  " << fieldName << nl;
    }


    template<>
    void Foam::rawSurfaceWriter::writeHeader
    (
        Ostream& os,
        const word& fieldName,
        const Field<vector>& values
    )
    {
        os  << values.size() << nl
            << "#  x  y  z"
            << "  " << fieldName << "_x"
            << "  " << fieldName << "_y"
            << "  " << fieldName << "_z"
            << nl;
    }


    template<>
    void Foam::rawSurfaceWriter::writeHeader
    (
        Ostream& os,
        const word& fieldName,
        const Field<sphericalTensor>& values
    )
    {
        os  << values.size() << nl
            << "#  ii  "
            << fieldName << "_ii" << nl;
    }


    template<>
    void Foam::rawSurfaceWriter::writeHeader
    (
        Ostream& os,
        const word& fieldName,
        const Field<symmTensor>& values
    )
    {
        os  << values.size() << nl
            << "#  xx  xy  xz  yy  yz";
        for (int i=0; i<6; ++i)
        {
            os  << "  " << fieldName << "_" << int(i);
        }
        os  << nl;
    }


    template<>
    void Foam::rawSurfaceWriter::writeHeader
    (
        Ostream& os,
        const word& fieldName,
        const Field<tensor>& values
    )
    {
        os  << values.size() << nl
            << "#  xx  xy  xz  yx  yy  yz  zx  zy  zz";
        for (int i=0; i<9; ++i)
        {
            os  << "  " << fieldName << "_" << int(i);
        }
        os  << nl;
    }
}
#endif


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceWriters::rawSeq::rawSeq
(
    const pointField& points,
    const faceList& faces
)
:
    writerBase(points, faces),
    writeCompression_(IOstream::UNCOMPRESSED)
{}


Foam::surfaceWriters::rawSeq::rawSeq
(
    const pointField& points,
    const faceList& faces,
    const dictionary& options
)
:
    writerBase(points, faces),
    writeCompression_(IOstream::UNCOMPRESSED)
{
    if (options.found("compression"))
    {
        writeCompression_ =
            IOstream::compressionEnum(options.get<word>("compression"));
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

#if 0
Foam::fileName Foam::rawSurfaceWriter::write
(
    const fileName& surfaceName,
    const meshedSurf& surf,
    const bool verbose
) const
{
    // geometry:  rootdir/time/surfaceName.raw

    const fileName outputFile
    (
        outputDirectory()/timeName()/surfaceName + ".raw"
    );

    if (verbose)
    {
        Info<< "Writing geometry to " << outputFile << nl;
    }

    const pointField& points = surf.points();
    const faceList&    faces = surf.faces();


    if (!isDir(outputFile.path()))
    {
        mkDir(outputFile.path());
    }


    OFstream os
    (
        outputFile,
        IOstream::ASCII,
        IOstream::currentVersion,
        writeCompression_
    );

    // Header
    os  << "# geometry NO_DATA " << faces.size() << nl
        << "#  x  y  z" << nl;

    // Write faces centres
    for (const face& f : faces)
    {
        writeLocation(os, f.centre(points));
        os << nl;
    }

    os  << nl;

    return outputFile;
}
#endif


// ************************************************************************* //
