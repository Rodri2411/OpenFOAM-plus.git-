/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2018 OpenCFD Ltd.
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
    Foam::functionObjects::ensightWrite

Group
    grpUtilitiesFunctionObjects

Description
    Writes fields in ensight format.

    Example of function object specification:
    \verbatim
    ensight
    {
        type            ensightWrite;
        libs            ("libutilityFunctionObjects.so");
        writeControl    writeTime;
        writeInterval   1;
        format          binary;

        overwrite       true;
        width           12;

        fields          (U p);
    }
    \endverbatim

    \heading Basic Usage
    \table
        Property    | Description                           | Required | Default
        type        | Type name: ensightWrite               | yes |
        fields      | Fields to output                      | yes |
        boundary    | Convert boundary fields               | no  | true
        internal    | Convert internal fields               | no  | true
        nodeValues  | Write values at nodes                 | no  | false
    \endtable

    \heading Ensight Output Options
    \table
        Property    | Description                           | Required | Default
        format      | ascii or binary format        | no | same as simulation
        width       | Mask width for \c data/XXXX           | no  | 8
        directory   | The output directory name     | no | postProcessing/NAME
        overwrite   | Remove existing directory             | no  | false
        consecutive | Consecutive output numbering          | no  | false
        nodeValues  | Write values at nodes                 | no  | false
    \endtable

    \heading Output Selection
    \table
        Property    | Description                           | Required | Default
        region      | Name for a single region              | no  | region0
        zone        | Limit to cell zone (name or regex)    | no  |
        zones       | Limit to cell zones (names, regexs)   | no  |
        faceZones   | Select faceZones to write             | no  |
        patches     | Limit to listed patches (wordRe list) | no  |
        geometry    | Dictionary of enclosing volumes       | no  | empty dict
        bounds      | Clip with a bounding box              | no  |
    \endtable

Note
    A combination of zone names, enclosing volumes and bounding box clipping
    can be used to restrict the region of interest. The selection of
    enclosing volumes is treated as an OR operation. These selections are
    combined in an AND operation with any cell zones specified.
    The cell selections are clipped according to the bounding box.
    Omitting the patches entry is the same as specifying the conversion of all
    patches.

    Note that if the \c patches entry is an empty list, this will select all
    patches and suppress writing the internalMesh.
    Consecutive output numbering can be used in conjunction with \c overwrite.

See also
    Foam::functionObjects::vtkWrite
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::timeControl

SourceFiles
    ensightWrite.C
    ensightWriteTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_ensightWrite_H
#define functionObjects_ensightWrite_H

#include "fvMeshFunctionObject.H"
#include "ensightCase.H"
#include "ensightMesh.H"

#include "interpolation.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMeshSubsetProxy.H"
#include "searchableSurfaces.H"
#include "IOobjectNames.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class dictionary;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class ensightWrite Declaration
\*---------------------------------------------------------------------------*/

class ensightWrite
:
    public fvMeshFunctionObject
{
    // Private data

        //- Ensight writer options
        ensightMesh::options writeOpts_;

        //- Ensight case options
        ensightCase::options caseOpts_;

        //- The output directory
        fileName outputDir_;

        //- Consecutive output numbering
        bool consecutive_;

        //- Track changes in mesh geometry
        enum polyMesh::readUpdateState meshState_;

        //- Requested names of cellZones to process
        wordRes selectZones_;

        //- Requested names of fields to process
        wordRes selectFields_;

        //- Request conversion volume(s) based on searchable surfaces
        autoPtr<searchableSurfaces> selectEnclosed_;

        //- Restrict to bounding box
        boundBox bounds_;

        //- Mesh subset handler
        fvMeshSubset meshSubset_;

        //- Ensight case handler
        autoPtr<ensightCase> ensCase_;

        //- Ensight mesh handler
        autoPtr<ensightMesh> ensMesh_;


    // Private Member Functions

        //- Ensight case handler
        ensightCase& ensCase()
        {
            return *ensCase_;
        }

        //- Ensight mesh handler
        ensightMesh& ensMesh()
        {
            return *ensMesh_;
        }


        //- Load searchable surfaces for enclosing volumes
        autoPtr<searchableSurfaces> loadSearchableSurfaces
        (
            const dictionary& dict,
            const word& dictName
        ) const;

        //- Update mesh subset according to zones, geometry, bounds
        bool updateSubset(fvMeshSubset& subsetter) const;

        //- Read information for selections
        bool readSelection(const dictionary& dict);

        //- Update meshes, subMeshes etc.
        bool update();


    // Write

        //- Write all volume fields
        label writeAllVolFields
        (
            const fvMeshSubset& proxy,
            const IOobjectNames& allFieldNames
        );

        //- Write selected volume fields.
        template<class Type>
        label writeVolFields
        (
            const fvMeshSubset& proxy,
            const IOobjectNames& allFieldNames
        );


        //- No copy construct
        ensightWrite(const ensightWrite&) = delete;

        //- No copy assignment
        void operator=(const ensightWrite&) = delete;


public:

    //- Runtime type information
    TypeName("ensightWrite");


    // Constructors

        //- Construct from runTime and dictionary.
        ensightWrite
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~ensightWrite() = default;


    // Member Functions

        //- Read the ensightWrite specification
        virtual bool read(const dictionary& dict);

        //- Do nothing
        virtual bool execute();

        //- Write fields, flush case file
        virtual bool write();

        //- Do nothing at the final time-loop
        virtual bool end();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh& mpm);

        //- Update for mesh point-motion
        virtual void movePoints(const polyMesh& mpm);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ensightWriteTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
