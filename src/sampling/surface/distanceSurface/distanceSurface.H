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

Class
    Foam::distanceSurface

Description
    A surface defined by distance to an input surface - using either
    an isoSurfaceCell or an isoSurface.

Usage
    Dictionary controls:
    \table
        Property | Description                             | Required | Default
        distance | distance from surface                   | yes      |
        signed   | apply sign for +ve/-ve distance         | yes      |
        cell     | use isoSurfaceCell algorithm            | no       | true
        regularise | point snapping                        | yes      |
        bounds   | limit with bounding box                 | no       |
        surfaceType | type of surface                      | yes      |
        surfaceName | name of surface in triSurface/       | no  | dict name
    \endtable


SourceFiles
    distanceSurface.C

\*---------------------------------------------------------------------------*/

#ifndef distanceSurface_H
#define distanceSurface_H

#include "sampledSurface.H"
#include "searchableSurface.H"
#include "isoSurfaceCell.H"
#include "isoSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class distanceSurface Declaration
\*---------------------------------------------------------------------------*/

class distanceSurface
{
    // Private data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Surface
        const autoPtr<searchableSurface> surfPtr_;

        //- Distance value
        const scalar distance_;

        //- Signed distance
        const bool signed_;

        //- Whether to use isoSurfaceCell or isoSurface
        const bool cell_;

        //- Whether to coarsen
        const bool regularise_;

        //- Optional bounding box to trim triangles against
        const boundBox bounds_;

        //- If restricted to zones, name of this zone or a regular expression
        keyType zoneKey_;

        //- Distance to cell centres
        autoPtr<volScalarField> cellDistancePtr_;

        //- Distance to points
        scalarField pointDistance_;

        //- Constructed iso surface
        autoPtr<isoSurfaceCell> isoSurfCellPtr_;

        //- Constructed iso surface
        autoPtr<isoSurface> isoSurfPtr_;


public:

    //- Runtime type information
    TypeName("distanceSurface");


    // Constructors

        //- Construct from dictionary
        distanceSurface
        (
            const word& defaultSurfaceName,
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Construct from components
        distanceSurface
        (
            const polyMesh& mesh,
            const bool interpolate,
            const word& surfaceType,
            const word& surfaceName,
            const scalar distance,
            const bool signedDistance,
            const bool cell,
            const bool regularise,
            const boundBox& bounds = boundBox::invertedBox
        );


    //- Destructor
    virtual ~distanceSurface() = default;


    // Member Functions

        //- Create/recreate the distance surface
        void createGeometry();

        //- The name of the underlying searchableSurface
        const word& surfaceName() const
        {
            return surfPtr_->name();
        }

        //- The distance to the underlying searchableSurface
        scalar distance() const
        {
            return distance_;
        }


        //- The underlying surface
        const meshedSurface& surface() const
        {
            if (cell_)
            {
                return *isoSurfCellPtr_;
            }
            return *isoSurfPtr_;
        }


        //- The underlying surface
        meshedSurface& surface()
        {
            if (cell_)
            {
                return *isoSurfCellPtr_;
            }
            return *isoSurfPtr_;
        }

        //- For each face, the original cell in mesh
        const labelList& meshCells() const
        {
            if (cell_)
            {
                return isoSurfCellPtr_->meshCells();
            }
            return isoSurfPtr_->meshCells();
        }

        //- For each face, the original cell in mesh
        labelList& meshCells()
        {
            if (cell_)
            {
                return isoSurfCellPtr_->meshCells();
            }
            return isoSurfPtr_->meshCells();
        }


    // Interpolate

        //- Interpolate volume field onto surface points
        template<class Type>
        tmp<Field<Type>> interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>& cellValues,
            const Field<Type>& pointValues
        ) const;


    // Output

        //- Print information
        void print(Ostream& os) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "distanceSurfaceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
