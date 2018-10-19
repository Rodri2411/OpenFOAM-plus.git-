/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2018 OpenCFD Ltd.
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

#include "sampledSurfaces.H"
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"
#include "IOmanip.H"
#include "interpolationCell.H"
#include "volPointInterpolation.H"
#include "PatchTools.H"
#include "mapPolyMesh.H"
#include "sampledTriSurfaceMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledSurfaces, 1);

    addToRunTimeSelectionTable
    (
        functionObject,
        sampledSurfaces,
        dictionary
    );
}

bool Foam::sampledSurfaces::verbose_ = false;
Foam::scalar Foam::sampledSurfaces::mergeTol_ = 1e-10;

const Foam::wordList Foam::sampledSurfaces::fieldTypeNames
({
    Foam::volScalarField::typeName,
    Foam::volVectorField::typeName,
    Foam::volSphericalTensorField::typeName,
    Foam::volSymmTensorField::typeName,
    Foam::volTensorField::typeName,

    Foam::volScalarField::Internal::typeName,
    Foam::volVectorField::Internal::typeName,
    Foam::volSphericalTensorField::Internal::typeName,
    Foam::volSymmTensorField::Internal::typeName,
    Foam::volTensorField::Internal::typeName,

    Foam::surfaceScalarField::typeName,
    Foam::surfaceVectorField::typeName,
    Foam::surfaceSphericalTensorField::typeName,
    Foam::surfaceSymmTensorField::typeName,
    Foam::surfaceTensorField::typeName
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sampledSurfaces::writeGeometry() const
{
    // Write to time directory under outputPath_
    // Skip surfaces without faces (eg, a failed cut-plane)

    for (sampledSurfaceProxy& proxy : proxies_)
    {
        if (proxy.hasGeometry())
        {
            formatter_->write(proxy.name(), proxy.geometry());

            proxy.handled();
        }
    }
}


void Foam::sampledSurfaces::writeOriginalIds()
{
    const word fieldName = "Ids";

    for (sampledSurfaceProxy& proxy : proxies_)
    {
        const sampledSurface& s = proxy.surface();

        if (s.hasFaceIds())
        {
            // Transcribe from label to scalar
            Field<scalar> ids
            (
                ListOps::create<scalar>
                (
                    s.originalIds(),
                    [](const label& val) -> scalar { return scalar(val); }
                )
            );

            writeSurface<scalar>(proxy, fieldName, ids);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledSurfaces::sampledSurfaces
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObjects::regionFunctionObject(name, runTime, dict),
    PtrList<sampledSurface>(),
    mesh_(refCast<const fvMesh>(obr_)),
    loadFromFiles_(false),
    outputPath_
    (
        time_.globalPath()/functionObject::outputPrefix/name
    ),
    fieldSelection_(),
    sampleFaceScheme_(),
    sampleNodeScheme_(),
    proxies_(),
    formatter_(nullptr)
{
    outputPath_.clean();  // Remove unneeded ".."

    read(dict);
}


Foam::sampledSurfaces::sampledSurfaces
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    functionObjects::regionFunctionObject(name, obr, dict),
    PtrList<sampledSurface>(),
    mesh_(refCast<const fvMesh>(obr)),
    loadFromFiles_(loadFromFiles),
    outputPath_
    (
        time_.globalPath()/functionObject::outputPrefix/name
    ),
    fieldSelection_(),
    sampleFaceScheme_(),
    sampleNodeScheme_(),
    proxies_(),
    formatter_(nullptr)
{
    outputPath_.clean();  // Remove unneeded ".."

    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sampledSurfaces::verbose(const bool verbosity)
{
    verbose_ = verbosity;
}


bool Foam::sampledSurfaces::execute()
{
    return true;
}


bool Foam::sampledSurfaces::write()
{
    if (empty())
    {
        return true;
    }

    // Finalize surfaces, merge points etc.
    update();

    // Update time values for output
    formatter_->setTime(time_.value(), time_.timeName());

    IOobjectNames objNames(classifyFields());

    // Write geometry first if required,
    // or when no fields would otherwise be written
    if (formatter_->separateGeometry() || objNames.empty())
    {
        writeGeometry();
    }

    const IOobjectList objects(obr_, obr_.time().timeName());

    sampleAndWrite<volScalarField>(objects);
    sampleAndWrite<volVectorField>(objects);
    sampleAndWrite<volSphericalTensorField>(objects);
    sampleAndWrite<volSymmTensorField>(objects);
    sampleAndWrite<volTensorField>(objects);

    sampleAndWrite<surfaceScalarField>(objects);
    sampleAndWrite<surfaceVectorField>(objects);
    sampleAndWrite<surfaceSphericalTensorField>(objects);
    sampleAndWrite<surfaceSymmTensorField>(objects);
    sampleAndWrite<surfaceTensorField>(objects);

    return true;
}


bool Foam::sampledSurfaces::read(const dictionary& dict)
{
    if (dict.found("surfaces"))
    {
        sampleFaceScheme_ = dict.lookupOrDefault<word>("sampleScheme", "cell");

        dict.readEntry("interpolationScheme", sampleNodeScheme_);
        dict.readEntry("fields", fieldSelection_);

        const word writeType(dict.get<word>("surfaceFormat"));

        // Define the surface formatter
        // Optionally defined extra controls for the output formats
        formatter_ = surfaceWriter::New
        (
            writeType,
            dict.subOrEmptyDict("formatOptions").subOrEmptyDict(writeType)
        );

        // Set the output directory
        formatter_->outputDirectory(outputPath_);

        // Set the output time value
        formatter_->setTime(time_.value(), time_.timeName());


        PtrList<sampledSurface>& surfs = surfaces();

        surfs = PtrList<sampledSurface>
        (
            dict.lookup("surfaces"),
            sampledSurface::iNew(mesh_)
        );

        proxies_.resize(surfs.size());

        forAll(proxies_, surfi)
        {
            proxies_.set
            (
                surfi,
                new sampledSurfaceProxy(surfs[surfi])
            );
        }

        // Ensure all surfaces and merge information are properly expired
        expire();

        if (surfs.size())
        {
            Info<< "Reading surface description:" << nl;
            for (const sampledSurface& s : surfs)
            {
                Info<< "    " << s.name() << nl;
            }
            Info<< endl;
        }
    }

    if (Pstream::master() && debug)
    {
        Pout<< "sample fields:" << fieldSelection_ << nl
            << "sample surfaces:" << nl << "(" << nl;

        for (const sampledSurface& s : surfaces())
        {
            Pout<< "  " << s << nl;
        }
        Pout<< ")" << endl;
    }

    return true;
}


void Foam::sampledSurfaces::updateMesh(const mapPolyMesh& mpm)
{
    if (&mpm.mesh() == &mesh_)
    {
        expire();
    }

    // pointMesh and interpolation will have been reset in mesh.update
}


void Foam::sampledSurfaces::movePoints(const polyMesh& mesh)
{
    if (&mesh == &mesh_)
    {
        expire();
    }
}


void Foam::sampledSurfaces::readUpdate(const polyMesh::readUpdateState state)
{
    if (state != polyMesh::UNCHANGED)
    {
        expire();
    }
}


bool Foam::sampledSurfaces::needsUpdate() const
{
    for (const sampledSurface& s : surfaces())
    {
        if (s.needsUpdate())
        {
            return true;
        }
    }

    return false;
}


bool Foam::sampledSurfaces::expire()
{
    label nChanged = 0;

    for (sampledSurfaceProxy& proxy : proxies_)
    {
        if (proxy.expire())
        {
            ++nChanged;
        }
    }

    // True if anything just expired
    return nChanged;
}


bool Foam::sampledSurfaces::update()
{
    if (!needsUpdate())
    {
        return false;
    }

    label nChanged = 0;


    // Dimension as fraction of mesh bounding box
    const scalar mergeDim = mergeTol_*mesh_.bounds().mag();

    if (Pstream::parRun())
    {
        DebugInfo
            << nl << "Merging all points within "
            << mergeDim << " metre" << nl;
    }

    for (sampledSurfaceProxy& proxy : proxies_)
    {
        if (proxy.update(mergeDim))
        {
            ++nChanged;
        }
    }

    return nChanged;
}


Foam::scalar Foam::sampledSurfaces::mergeTol()
{
    return mergeTol_;
}


Foam::scalar Foam::sampledSurfaces::mergeTol(const scalar tol)
{
    const scalar prev(mergeTol_);
    mergeTol_ = tol;
    return prev;
}


// ************************************************************************* //
