/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "coordinateSystem.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::coordinateSystem> Foam::coordinateSystem::New
(
    const objectRegistry& obr,
    const dictionary& dict
)
{
    const dictionary& coordDict = dict.subDict(typeName_());
    const word modelType(coordDict.get<word>("type"));

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "Unknown coordinateSystem type "
            << modelType << nl << nl
            << "Valid types:  "
            << flatOutput(dictionaryConstructorTablePtr_->sortedToc())
            << exit(FatalIOError);
    }

    return autoPtr<coordinateSystem>(cstrIter()(obr, coordDict));
}


Foam::autoPtr<Foam::coordinateSystem> Foam::coordinateSystem::New
(
    const dictionary& dict
)
{
    const dictionary& coordDict = dict.subDict(typeName_());

    return autoPtr<coordinateSystem>::New(coordDict);
}


Foam::autoPtr<Foam::coordinateSystem> Foam::coordinateSystem::New
(
    Istream& is
)
{
    const word name(is);
    const dictionary dict(is);

    return autoPtr<coordinateSystem>::New(name, dict);
}


// ************************************************************************* //
