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
    Foam::sampledSurfaceProxy

Description
    Management of sampled surfaces with merging

SourceFiles
    sampledSurfaceProxy.C

\*---------------------------------------------------------------------------*/

#ifndef sampledSurfaceProxy_H
#define sampledSurfaceProxy_H

#include <functional>
#include "mergedSurf.H"
#include "sampledNull.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class sampledSurfaceProxy Declaration
\*---------------------------------------------------------------------------*/

class sampledSurfaceProxy
{
public:

    //- Enumeration defining topology state.
    using updateState = sampledSurface::updateState;


private:

    // Static data members

        //- A placeholder surface
        static sampledNull dummy_;


    // Private data

        //- Reference to the original sampledSurface
        std::reference_wrapper<sampledSurface> surf_;

        //- Merging information
        mergedSurf merged_;

        //- The topology state
        updateState state_;


public:

    // Constructors

        //- Construct null
        sampledSurfaceProxy();

        //- Construct for the given surface, with default merge distance
        explicit sampledSurfaceProxy(sampledSurface& s);

        //- Construct for the given surface, with merge distance
        sampledSurfaceProxy(sampledSurface& s, scalar mergeDim);


    //- Destructor
    virtual ~sampledSurfaceProxy() = default;


    // Member Functions

    // Access

        //- The surface name
        const word& name() const
        {
            return surface().name();
        }

        //- Reference to the original sampledSurface
        const sampledSurface& surface() const;

        //- Reference to the original sampledSurface
        sampledSurface& surface();


        //- True if the surface or merged surface (master only) has faces
        bool hasGeometry() const
        {
            return
            (
                Pstream::parRun()
              ? (Pstream::master() && !merged_.faces().empty())
              : !surface().faces().empty()
            );
        }


        //- The effective surface geometry
        const meshedSurf& geometry() const
        {
            if (Pstream::parRun())
            {
                return merged_;
            }

            return surf_.get();
        }


        //- Return the current updateState
        updateState state() const
        {
            return state_;
        }


    // Edit

        //- Change reference to the original sampledSurface
        void surface(sampledSurface& s);

        //- Merge points
        void merge(scalar mergeDim);

        //- Mark the surface as needing an update.
        //  May also free up unneeded data.
        //  Return false if surface was already marked as expired.
        bool expire();

        //- Update the surface as required and merge the resulting points.
        //  Do nothing (and return false) if no update was required
        bool update(scalar mergeDim);

        //- Any geometry changes pending?
        bool staged() const;

        //- Handled the geometry changes (ie, in a writer)
        bool handled();


    // Other

        //- Renumber (point data) to correspond to merged points
        template<class Type>
        bool renumberPointData(List<Type>& values)
        {
            const bool ok = (merged_.pointsMap().size() == values.size());

            if (ok)
            {
                inplaceReorder(merged_.pointsMap(), values);
                values.resize(merged_.points().size());
            }

            return ok;
        }
};


// Global Operators

//- Ostream operator
Ostream& operator<<(Ostream& os, const sampledSurfaceProxy& proxy);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
