/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::surfZoneIOList

Description
    IOobject for a surfZoneList

SourceFiles
    surfZoneIOList.C

\*---------------------------------------------------------------------------*/

#ifndef surfZoneIOList_H
#define surfZoneIOList_H

#include "surfZoneList.H"
#include "regIOobject.H"
#include "faceList.H"
#include "className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations

class surfZoneIOList;

Ostream& operator<<(Ostream&, const surfZoneIOList&);


/*---------------------------------------------------------------------------*\
                      Class surfZoneIOList Declaration
\*---------------------------------------------------------------------------*/

class surfZoneIOList
:
    public regIOobject,
    public surfZoneList
{
public:

    //- Runtime type information
    TypeName("surfZoneList");


    // Constructors

        //- Construct from IOobject
        explicit surfZoneIOList(const IOobject& io);

        //- Construct from IOobject and surfZoneList content
        surfZoneIOList(const IOobject& io, const UList<surfZone>& content);

        //- Construct from IOobject and transferring surfZoneList content
        surfZoneIOList(const IOobject& io, surfZoneList&& content);


    //- Destructor
    virtual ~surfZoneIOList() = default;


    // Member Functions

        //- The writeData function for regIOobject
        bool writeData(Ostream& os) const;


    // Member Operators

        //- Copy assignment of entries
        void operator=(const surfZoneIOList& rhs);

        //- Copy or move assignment of entries
        using surfZoneList::operator=;


    // IOstream Operators

        friend Ostream& operator<<(Ostream& os, const surfZoneIOList& list);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
