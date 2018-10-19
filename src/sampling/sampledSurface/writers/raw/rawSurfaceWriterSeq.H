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

Class
    Foam::surfaceWriters:raw_serial

Description
    A surfaceWriter for raw output.

    \heading Output file locations

    The \c rootdir normally corresponds to something like
    \c postProcessing/\<name\>

    \subheading Geometry
    \verbatim
    rootdir
    `-- timeName
        `-- surfaceName.{raw}
    \endverbatim

    \subheading Fields
    \verbatim
    rootdir
    `-- timeName
        |-- <field0>_surfaceName.{raw}
        `-- <field1>_surfaceName.{raw}
    \endverbatim

SourceFiles
    rawSurfaceWriter.C

\*---------------------------------------------------------------------------*/

#ifndef rawSurfaceWriterSeq_H
#define rawSurfaceWriterSeq_H

#include "surfaceWriterBase.H"
#include "surfaceWritersFactory.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceWriters
{

/*---------------------------------------------------------------------------*\
                   Class surfaceWriters::rawSeq Declaration
\*---------------------------------------------------------------------------*/

class rawSeq
:
    public writerBase
{
    // Private data

        //- Output compression, defaults to uncompressed
        IOstream::compressionType writeCompression_;


    // Private Member Functions

        //- Emit x,y,z
        static inline void writeLocation(Ostream& os, const point& p)
        {
            os << p.x() << ' ' << p.y() << ' ' << p.z();
        }

        //- Emit each component
        template<class Type>
        static inline void writeData(Ostream& os, const Type& val)
        {
            const direction ncmpt = pTraits<Type>::nComponents;
            for (direction cmpt=0; cmpt < ncmpt; ++cmpt)
            {
                os  << ' ' << component(val, cmpt);
            }
            os  << nl;
        }


        template<class Type>
        static void writeHeader
        (
            Ostream& os,
            const word& fieldName,
            const Field<Type>& values
        );

//         //- Templated write operation
//         template<class Type>
//         fileName writeTemplate
//         (
//             const fileName& surfaceName,    //!< Name of surface
//             const meshedSurf& surf,         //!< Surface geometry
//             const word& fieldName,          //!< Name of field
//             const Field<Type>& values,      //!< Field values to write
//             const bool isNodeValues,        //!< Values are per-vertex
//             const bool verbose = false      //!< Additional verbosity
//         ) const;


public:

    // Public classes
    class factory;

    //- Runtime type information
    TypeName("raw");


    // Constructors

        //- Construct from components (default format INLINE_BASE64),
        //- and open the file for writing.
        //  The file name is with/without an extension.
        rawSeq
        (
            const pointField& points,
            const faceList& faces
        );

        //- Construct from components and open the file for writing.
        //  The file name is with/without an extension.
        rawSeq
        (
            const pointField& points,
            const faceList& faces,
            const dictionary& options
        );


    //- Destructor
    virtual ~rawSeq() = default;


    // Member Functions


        //- Write single surface geometry to file.
        virtual bool writeGeometry() { return true; };  // override

        virtual bool beginCellData() { return true; };  // override
        virtual bool beginPointData() { return true; };  // override

        virtual bool write
        (
            const word& fieldName, const UList<scalar>& field
        ) { return true; };  // override
        virtual bool write
        (
            const word& fieldName, const UList<vector>& field
        ) { return true; };  // override
        virtual bool write
        (
            const word& fieldName, const UList<sphericalTensor>& field
        ) { return true; };  // override
        virtual bool write
        (
            const word& fieldName, const UList<symmTensor>& field
        ) { return true; };  // override

        virtual bool write
        (
            const word& fieldName, const UList<tensor>& field
        ) { return true; };  // override

};


/*---------------------------------------------------------------------------*\
                   Class surfaceWriters::rawSeq Declaration
\*---------------------------------------------------------------------------*/

namespace factories
{

struct rawSeq
:
    public surfaceWriters::factory
{

    //- Construct null
    rawSeq()
    :
        surfaceWriters::factory()
    {}

    //- Construct with options
    rawSeq(const dictionary& options)
    :
        surfaceWriters::factory(options)
    {}


    //- Create a new writer
    virtual autoPtr<surfaceWriters::writerBase> newWriter
    (
        const pointField& points,
        const faceList& faces
    ) const
    {
        return
            autoPtr<surfaceWriters::writerBase>::NewFrom
            <surfaceWriters::rawSeq>(points, faces, options());
    }

    //- Create a new writer
    virtual autoPtr<surfaceWriters::writerBase> newWriter
    (
        const pointField& points,
        const faceList& faces,
        const dictionary& options
    ) const
    {
        return
            autoPtr<surfaceWriters::writerBase>::NewFrom
            <surfaceWriters::rawSeq>(points, faces, options);
    }
};

} // End namespace surfaceWriters


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceWriters
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
//     #include "rawSurfaceWriterSeqTemplates.C"
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
