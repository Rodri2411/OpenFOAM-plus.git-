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

#include "foamSurfaceWriter.H"
#include "OFstream.H"
#include "makeSurfaceWriterMethods.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    makeSurfaceWriterType(foamSurfaceWriter);

    // Field writing methods
    defineSurfaceWriterWriteFields(foamSurfaceWriter);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Field writing implementation
#include "foamSurfaceWriterImpl.C"

// Field writing methods
defineSurfaceWriterWriteFields(Foam::foamSurfaceWriter);


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::foamSurfaceWriter::write
(
    const fileName& surfaceName,
    const meshedSurf& surf,
    const bool verbose
) const
{
    // Output:
    // - rootdir/time/surfaceName/{points,faces}

    const fileName base(outputDirectory() / timeName() / surfaceName);

    if (!isDir(base))
    {
        mkDir(base);
    }

    if (verbose)
    {
        Info<< "Writing geometry to " << base << endl;
    }

    const pointField& points = surf.points();
    const faceList&    faces = surf.faces();

    // Points
    {
        OFstream(base/"points")() << points;
    }

    // Faces
    {
        OFstream(base/"faces")() << faces;
    }

    // Face centers.
    // Not really necessary but very handy when reusing as inputs,
    // e.g. for timeVaryingMapped BC.
    {
        pointField faceCentres(faces.size(), Zero);

        forAll(faces, facei)
        {
            faceCentres[facei] = faces[facei].centre(points);
        }

        OFstream(base/"faceCentres")() << faceCentres;
    }

    return base;
}


// ************************************************************************* //
