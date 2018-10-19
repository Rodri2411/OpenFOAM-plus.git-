/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017-2018 OpenCFD Ltd.
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
#include "IOobjectList.H"
#include "stringListOps.H"
#include "UIndirectList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::IOobjectNames Foam::sampledSurfaces::classifyFields() const
{
    wordList allFields;    // For warnings

    IOobjectNames objNames;

    if (loadFromFiles_)
    {
        // Check files for a particular time
        IOobjectList objects(obr_, obr_.time().timeName());

        allFields = objects.sortedNames();
        objNames = IOobjectNames(objects, fieldSelection_);
    }
    else
    {
        // Check currently available fields
        allFields = obr_.sortedNames();
        objNames  = IOobjectNames(obr_, fieldSelection_);
    }

    DynamicList<word> missed(fieldSelection_.size());

    // Detect missing fields
    for (const wordRe& pattern : fieldSelection_)
    {
        if (pattern.isLiteral() && !allFields.found(pattern))
        {
            missed.append(pattern);
        }
    }

    if (missed.size())
    {
        WarningInFunction
            << nl
            << "Cannot find specified "
            << (loadFromFiles_ ? "field file(s)" : "registered field(s)")
            << nl << "    " << flatOutput(missed) << nl << endl;
    }

    objNames.retainClasses(fieldTypeNames);

    return objNames;
}


// ************************************************************************* //
