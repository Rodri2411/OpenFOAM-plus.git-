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

#include "surfaceWritersFactory.H"
#include "HashTable.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceWriters
{
    defineTypeNameAndDebug(factory, 0);
    defineRunTimeSelectionTable(factory, word);
    defineRunTimeSelectionTable(factory, wordDict);
//     addNamedToRunTimeSelectionTable
//     (
//         factory,
//         factory,
//         word,
//         null
//     );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceWriters::factory::factory()
{}


Foam::surfaceWriters::factory::factory(const dictionary& options)
:
    options_(options)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::surfaceWriters::factory>
Foam::surfaceWriters::factory::New(const word& writerType)
{
    auto cstrIter = wordConstructorTablePtr_->cfind(writerType);

//     if (!cstrIter.found())
//     {
//        if (MeshedSurfaceProxy<face>::canWriteType(writerType))
//        {
//            // generally unknown, but can be written via MeshedSurfaceProxy
//            // use 'proxy' handler instead
//            return autoPtr<surfaceWriter>(new proxySurfaceWriter(writerType));
//        }
//
    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown writer type \"" << writerType << "\"\n\n"
            << "Valid write types : "
            << wordConstructorTablePtr_->sortedToc() << nl
            << exit(FatalError);
    }

    return autoPtr<surfaceWriters::factory>(cstrIter()());
}


Foam::autoPtr<Foam::surfaceWriters::factory>
Foam::surfaceWriters::factory::New
(
    const word& writerType,
    const dictionary& writerOptions
)
{
    // Find constructors with dictionary options
    auto cstrIter = wordDictConstructorTablePtr_->cfind(writerType);

    if (!cstrIter.found())
    {
        // Revert to versions without options
        return surfaceWriters::factory::New(writerType);
    }

    return autoPtr<surfaceWriters::factory>(cstrIter()(writerOptions));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::dictionary& Foam::surfaceWriters::factory::options() const
{
    return options_;
}


// ************************************************************************* //
