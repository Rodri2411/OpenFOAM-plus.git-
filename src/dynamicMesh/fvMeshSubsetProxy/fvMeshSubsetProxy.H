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
    Foam::fvMeshSubsetProxy

Description
    Simple proxy for holding a mesh, or mesh-subset.
    The subMeshes are currently limited to cellSet or cellZone definitions.

SourceFiles
    fvMeshSubsetProxy.C
    fvMeshSubsetProxyTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshSubsetProxy_H
#define fvMeshSubsetProxy_H

#include "wordRes.H"
#include "fvMeshSubset.H"
#include "zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class Time;

/*---------------------------------------------------------------------------*\
                      Class fvMeshSubsetProxy Declaration
\*---------------------------------------------------------------------------*/

class fvMeshSubsetProxy
{
public:

        //- Internal bookkeeping for subset type
        enum subsetType
        {
            NONE = 0,       //<! No subset
            SET  = 0x1,     //<! Subset with a cellSet
            ZONE = 0x2,     //<! Subset with a cellZone
            ZONES = 0x4     //<! Subset with multiple cellZones
        };


private:

    // Private data

        //- Reference to mesh
        fvMesh& baseMesh_;

        //- Subsetting engine
        fvMeshSubset subsetter_;

        //- Patch ID for exposed internal faces
        label exposedPatchId_;

        //- Subsetting type. From subsetType enum and bounds flag.
        unsigned type_;

        //- Name of the cellSet/cellZone (or empty)
        const word name_;

        //- Selection for multiple cell zones
        wordRes zones_;

        //- Optional bounding box
        boundBox bounds_;


    // Private Member Functions

        //- No copy construct
        fvMeshSubsetProxy(const fvMeshSubsetProxy&) = delete;

        //- No copy assignment
        void operator=(const fvMeshSubsetProxy&) = delete;


public:

    // Constructors

        //- Construct a pass-through proxy. No correct() invoked or required.
        explicit fvMeshSubsetProxy(fvMesh& baseMesh);

        //- Construct from mesh with bounding box limiter.
        //  Invokes correct() if the bounding is non-empty.
        fvMeshSubsetProxy(fvMesh& baseMesh, const boundBox& bounds);

        //- Construct from components
        //  Optionally defer invoking correct(), which can be useful
        //  when adding a bounding box
        fvMeshSubsetProxy
        (
            fvMesh& baseMesh,
            const subsetType,
            const word& name,
            const boundBox& bounds = boundBox::invertedBox,
            label exposedPatchId = -1
        );

        //- Construct from components
        fvMeshSubsetProxy
        (
            fvMesh& baseMesh,
            const wordRes& zones,
            const word& name,
            const boundBox& bounds = boundBox::invertedBox,
            label exposedPatchId = -1
        );


    // Member Functions

    // Access

        //- The entire base mesh
        inline const fvMesh& baseMesh() const
        {
            return baseMesh_;
        }

        //- The mesh subsetter
        inline const fvMeshSubset& subsetter() const
        {
            return subsetter_;
        }

        //- Check if a sub-mesh is being used
        inline bool useSubMesh() const
        {
            return type_ != NONE;
        }

        //- Access either base-mesh or sub-mesh
        inline const fvMesh& mesh() const
        {
            if (useSubMesh())
            {
                return subsetter_.subMesh();
            }

            return baseMesh_;
        }

        //- The associated (set or zone) name if any.
        inline const word& name() const
        {
            return name_;
        }


    // Edit

        //- Define/change bounding limits.
        void bounds(const boundBox& bounds);

        //- Update of mesh subset
        void correct(bool verbose = false);

        //- Read mesh. Correct on topo-change
        polyMesh::readUpdateState readUpdate();


    // Fields

        //- Construct volField (with zeroGradient) from an internal field
        template<class Type>
        static tmp<GeometricField<Type, fvPatchField, volMesh>>
        zeroGradientField
        (
            const DimensionedField<Type, volMesh>& df
        );

        //- Convert an internal field to a volume field (with zeroGradient)
        template<class Type>
        static tmp<GeometricField<Type, fvPatchField, volMesh>>
        interpolateInternal
        (
            const fvMeshSubset& subsetter,
            const DimensionedField<Type, volMesh>& df
        );

        //- Convert an internal field to a volume field (with zeroGradient)
        //  Currently no proper memory reuse
        template<class Type>
        static tmp<GeometricField<Type, fvPatchField, volMesh>>
        interpolateInternal
        (
            const fvMeshSubset& subsetter,
            const tmp<DimensionedField<Type, volMesh>>& tdf
        );

        //- Wrapper for field or the subsetted field.
        //  Pass through or forward to fvMeshSubset::interpolate()
        template<class GeoField>
        static tmp<GeoField> interpolate
        (
            const fvMeshSubset& subsetter,
            const GeoField& fld
        );

        //- Wrapper for field or the subsetted field.
        //  Pass through or forward to fvMeshSubset::interpolate()
        template<class GeoField>
        static tmp<GeoField> interpolate
        (
            const fvMeshSubset& subsetter,
            const tmp<GeoField>& fld
        );


    // Fields

        //- Convert an internal field to a volume field (with zeroGradient)
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh>>
        interpolateInternal
        (
            const DimensionedField<Type, volMesh>& df
        ) const;

        //- Convert an internal field to a volume field (with zeroGradient)
        //  Currently no proper memory reuse
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh>>
        interpolateInternal
        (
            const tmp<DimensionedField<Type, volMesh>>& tdf
        ) const;


        //- Wrapper for field or the subsetted field.
        //  Pass through or forward to fvMeshSubset::interpolate()
        template<class GeoField>
        tmp<GeoField> interpolate(const GeoField& fld) const;

        //- Wrapper for field or the subsetted field.
        //  Pass through or forward to fvMeshSubset::interpolate()
        template<class GeoField>
        tmp<GeoField> interpolate(const tmp<GeoField>& fld) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvMeshSubsetProxyTemplates.C"
#endif

#endif

// ************************************************************************* //
