/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "distanceSurface.H"
#include "dictionary.H"
#include "volFields.H"
#include "volPointInterpolation.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(distanceSurface, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distanceSurface::distanceSurface
(
    const word& defaultSurfaceName,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    surfPtr_
    (
        searchableSurface::New
        (
            dict.lookup("surfaceType"),
            IOobject
            (
                dict.lookupOrDefault("surfaceName", defaultSurfaceName),
                mesh.time().constant(), // directory
                "triSurface",           // instance
                mesh.time(),            // registry
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            dict
        )
    ),
    distance_(readScalar(dict.lookup("distance"))),
    signed_(dict.get<bool>("signed")),
    cell_(dict.lookupOrDefault("cell", true)),
    regularise_(dict.lookupOrDefault("regularise", true)),
    bounds_(dict.lookupOrDefault("bounds", boundBox::invertedBox)),
    zoneKey_(keyType::null),
    isoSurfCellPtr_(nullptr),
    isoSurfPtr_(nullptr)
{}


Foam::distanceSurface::distanceSurface
(
    const polyMesh& mesh,
    const bool interpolate,
    const word& surfaceType,
    const word& surfaceName,
    const scalar distance,
    const bool signedDistance,
    const bool cell,
    const bool regularise,
    const boundBox& bounds
)
:
    mesh_(mesh),
    surfPtr_
    (
        searchableSurface::New
        (
            surfaceType,
            IOobject
            (
                surfaceName,            // name
                mesh.time().constant(), // directory
                "triSurface",           // instance
                mesh.time(),            // registry
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            dictionary()
        )
    ),
    distance_(distance),
    signed_(signedDistance),
    cell_(cell),
    regularise_(regularise),
    bounds_(bounds),
    zoneKey_(keyType::null),
    isoSurfCellPtr_(nullptr),
    isoSurfPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::distanceSurface::createGeometry()
{
    if (debug)
    {
        Pout<< "distanceSurface::createGeometry updating geometry." << endl;
    }

    // Clear any stored topologies
    isoSurfCellPtr_.clear();
    isoSurfPtr_.clear();

    const fvMesh& fvm = static_cast<const fvMesh&>(mesh_);

    // Distance to cell centres
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    cellDistancePtr_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "cellDistance",
                fvm.time().timeName(),
                fvm.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            fvm,
            dimensionedScalar(dimLength, Zero)
        )
    );
    volScalarField& cellDistance = *cellDistancePtr_;

    // Internal field
    {
        const pointField& cc = fvm.C();
        scalarField& fld = cellDistance.primitiveFieldRef();

        List<pointIndexHit> nearest;
        surfPtr_().findNearest
        (
            cc,
            scalarField(cc.size(), GREAT),
            nearest
        );

        if (signed_)
        {
            vectorField norms;
            surfPtr_().getNormal(nearest, norms);

            forAll(norms, i)
            {
                const point diff(cc[i] - nearest[i].hitPoint());

                fld[i] = sign(diff & norms[i]) * Foam::mag(diff);
            }
        }
        else
        {
            forAll(nearest, i)
            {
                fld[i] = Foam::mag(cc[i] - nearest[i].hitPoint());
            }
        }
    }

    volScalarField::Boundary& cellDistanceBf =
        cellDistance.boundaryFieldRef();

    // Patch fields
    {
        forAll(fvm.C().boundaryField(), patchi)
        {
            const pointField& cc = fvm.C().boundaryField()[patchi];
            fvPatchScalarField& fld = cellDistanceBf[patchi];

            List<pointIndexHit> nearest;
            surfPtr_().findNearest
            (
                cc,
                scalarField(cc.size(), GREAT),
                nearest
            );

            if (signed_)
            {
                vectorField norms;
                surfPtr_().getNormal(nearest, norms);

                forAll(norms, i)
                {
                    const point diff(cc[i] - nearest[i].hitPoint());

                    fld[i] = sign(diff & norms[i]) * Foam::mag(diff);
                }
            }
            else
            {
                forAll(nearest, i)
                {
                    fld[i] = Foam::mag(cc[i] - nearest[i].hitPoint());
                }
            }
        }
    }


    // On processor patches the mesh.C() will already be the cell centre
    // on the opposite side so no need to swap cellDistance.


    // Distance to points
    pointDistance_.setSize(fvm.nPoints());
    {
        const pointField& pts = fvm.points();

        List<pointIndexHit> nearest;
        surfPtr_().findNearest
        (
            pts,
            scalarField(pts.size(), GREAT),
            nearest
        );

        if (signed_)
        {
            vectorField norms;
            surfPtr_().getNormal(nearest, norms);

            forAll(norms, i)
            {
                const point diff(pts[i] - nearest[i].hitPoint());

                pointDistance_[i] = sign(diff & norms[i]) * Foam::mag(diff);
            }
        }
        else
        {
            forAll(nearest, i)
            {
                pointDistance_[i] = Foam::mag(pts[i] - nearest[i].hitPoint());
            }
        }
    }


    if (debug)
    {
        Pout<< "Writing cell distance:" << cellDistance.objectPath() << endl;
        cellDistance.write();
        pointScalarField pDist
        (
            IOobject
            (
                "pointDistance",
                fvm.time().timeName(),
                fvm.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            pointMesh::New(fvm),
            dimensionedScalar(dimLength, Zero)
        );
        pDist.primitiveFieldRef() = pointDistance_;

        Pout<< "Writing point distance:" << pDist.objectPath() << endl;
        pDist.write();
    }


    // Direct from cell field and point field.
    if (cell_)
    {
        isoSurfCellPtr_.reset
        (
            new isoSurfaceCell
            (
                fvm,
                cellDistance,
                pointDistance_,
                distance_,
                regularise_,
                bounds_
            )
        );
    }
    else
    {
        isoSurfPtr_.reset
        (
            new isoSurface
            (
                cellDistance,
                pointDistance_,
                distance_,
                regularise_,
                bounds_
            )
        );
    }

    if (debug)
    {
        print(Pout);
        Pout<< endl;
    }
}


void Foam::distanceSurface::print(Ostream& os) const
{
    os  << "  surface:" << surfaceName()
        << "  distance:" << distance()
        << "  faces:" << surface().surfFaces().size()
        << "  points:" << surface().points().size();
}


// ************************************************************************* //
