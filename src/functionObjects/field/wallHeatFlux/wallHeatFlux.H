/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::functionObjects::wallHeatFlux

Group
    grpForcesFunctionObjects

Description
    Calculates and write the heat-flux at wall patches as the
    volScalarField field 'wallHeatFlux'.

    All wall patches are included by default; to restrict the calculation to
    certain patches, use the optional 'patches' entry.

Usage
    Example of function object specification:
    \verbatim
    wallHeatFlux1
    {
        type        wallHeatFlux;
        libs        ("libfieldFunctionObjects.so");
        ...
        region      fluid;
        patches     (".*Wall");
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property | Description                | Required   | Default value
        type     | type name: wallHeatFlux    | yes        |
        patches  | list of patches to process | no         | all wall patches
        qr       | name of radiative heat flux field | no  | qr
        region   | region to be evaluated     | no         | default region
    \endtable

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile
    Foam::functionObjects::timeControl

SourceFiles
    wallHeatFlux.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_wallHeatFlux_H
#define functionObjects_wallHeatFlux_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldsFwd.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class wallHeatFlux Declaration
\*---------------------------------------------------------------------------*/

class wallHeatFlux
:
    public fvMeshFunctionObject,
    public writeFile
{

protected:

    // Protected data

        //- Optional list of patches to process
        labelHashSet patchSet_;

        //- Name of radiative heat flux name, default = qr
        word qrName_;


    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;

        //- Calculate the heat-flux
        void calcHeatFlux
        (
            const volScalarField& alpha,
            const volScalarField& he,
            volScalarField& wallHeatFlux
        );


private:

    // Private member functions

        //- No copy construct
        wallHeatFlux(const wallHeatFlux&) = delete;

        //- No copy assignment
        void operator=(const wallHeatFlux&) = delete;


public:

    //- Runtime type information
    TypeName("wallHeatFlux");


    // Constructors

        //- Construct from Time and dictionary
        wallHeatFlux
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );


    //- Destructor
    virtual ~wallHeatFlux();


    // Member Functions

        //- Read the wallHeatFlux data
        virtual bool read(const dictionary&);

        //- Calculate the wall heat-flux
        virtual bool execute();

        //- Write the wall heat-flux
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
