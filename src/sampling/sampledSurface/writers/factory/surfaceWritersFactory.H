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
    Foam::surfaceWriter

Description
    Base class for surface writers

SourceFiles
    surfaceWriter.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceWritersFactory_H
#define surfaceWritersFactory_H

#include "typeInfo.H"
#include "autoPtr.H"
#include "fileName.H"
#include "instant.H"
#include "MeshedSurface.H"
#include "surfaceWriter.H"
#include "surfaceWriterBase.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace surfaceWriters
{

/*---------------------------------------------------------------------------*\
                        Class surfaceWriter Declaration
\*---------------------------------------------------------------------------*/

class factory
{
    dictionary options_;

public:

    //- Runtime type information
    TypeName("factory");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            factory,
            word,
            (),
            ()
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            factory,
            wordDict,
            (
                const dictionary& writerOptions
            ),
            (writerOptions)
        );


    // Selectors

        //- Return a reference to the selected surfaceWriter
        static autoPtr<factory> New
        (
            const word& writerType
        );

        //- Return a reference to the selected surfaceWriter
        //  Select with extra write option
        static autoPtr<factory> New
        (
            const word& writerType,
            const dictionary& options
        );


    // Constructors

        //- Construct null
        factory();

        //- Construct with options
        factory(const dictionary& options);


    //- Destructor
    virtual ~factory() = default;


    // Member Functions

        //- The configuration options
        const dictionary& options() const;

        //- Write single surface geometry to file.
        virtual autoPtr<surfaceWriters::writerBase> newWriter
        (
            const pointField& points,
            const faceList& faces
        ) const = 0;

//        //- Write single surface geometry to file.
//        virtual autoPtr<surfaceWriters::writerBase> newWriter
//        (
//            const pointField& points,
//            const faceList& faces,
//            const dictionary& dict
//        ) const
//        {
//            return this->newWriter(points, faces);
//        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceWriters
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
