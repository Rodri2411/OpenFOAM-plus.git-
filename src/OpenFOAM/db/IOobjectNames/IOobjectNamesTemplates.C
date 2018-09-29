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

\*---------------------------------------------------------------------------*/

#include "IOobjectNames.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class UnaryPredicate>
Foam::label Foam::IOobjectNames::filter
(
    const UnaryPredicate& pred,
    bool pruning
)
{
    label changed = 0;

    forAllIters(classes_, iter)
    {
        changed += iter.object().filterKeys(pred, pruning);
    }

    return changed;
}


template<class ClassType>
Foam::wordList Foam::IOobjectNames::sortedNames() const
{
    auto iter = classes_.cfind(ClassType::typeName);

    if (iter.found())
    {
        return (*iter).sortedToc();
    }

    return wordList();
}


// ************************************************************************* //
