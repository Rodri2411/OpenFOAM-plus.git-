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
#include "IOobjectList.H"
#include "objectRegistry.H"
#include "HashOps.H"
#include "ListOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Container>
Foam::label Foam::IOobjectNames::countImpl
(
    const HashTable<wordHashSet>& classes,
    const Container& whichTypes
)
{
    label total = 0;

    for (const word& clsName : whichTypes)
    {
        auto iter = classes.cfind(clsName);

        if (iter.found())
        {
            total += (*iter).size();
        }
    }

    return total;
}


template<class Container>
void Foam::IOobjectNames::retainClassesImpl(const Container& whichTypes)
{
    wordHashSet others(classes_.toc());

    others.erase(whichTypes);
    classes_.erase(others);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::label Foam::IOobjectNames::count
(
    const objectRegistry& objects,
    const UList<word>& knownTypes,
    const wordRes& selection,
    const bool filter
)
{
    HashTable<wordHashSet> classes =
    (
        filter
      ? objects.classes(selection)
      : objects.classes()
    );

    return countImpl(classes, knownTypes);
}


Foam::label Foam::IOobjectNames::count
(
    const objectRegistry& objects,
    const wordHashSet& knownTypes,
    const wordRes& selection,
    const bool filter
)
{
    HashTable<wordHashSet> classes =
    (
        filter
      ? objects.classes(selection)
      : objects.classes()
    );

    return countImpl(classes, knownTypes);
}


Foam::label Foam::IOobjectNames::count
(
    const IOobjectList& objects,
    const UList<word>& knownTypes,
    const wordRes& selection,
    const bool filter
)
{
    HashTable<wordHashSet> classes =
    (
        filter
      ? objects.classes(selection)
      : objects.classes()
    );

    return countImpl(classes, knownTypes);
}


Foam::label Foam::IOobjectNames::count
(
    const IOobjectList& objects,
    const wordHashSet& knownTypes,
    const wordRes& selection,
    const bool filter
)
{
    HashTable<wordHashSet> classes =
    (
        filter
      ? objects.classes(selection)
      : objects.classes()
    );

    return countImpl(classes, knownTypes);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IOobjectNames::IOobjectNames()
:
    classes_()
{}


Foam::IOobjectNames::IOobjectNames(const objectRegistry& objects)
:
    classes_(objects.classes())
{}


Foam::IOobjectNames::IOobjectNames(const IOobjectList& objects)
:
    classes_(objects.classes())
{}


Foam::IOobjectNames::IOobjectNames
(
    const objectRegistry& objects,
    const wordRes& selection,
    const bool filter
)
:
    classes_
    (
        filter
      ? objects.classes(selection)
      : objects.classes()
    )
{}


Foam::IOobjectNames::IOobjectNames
(
    const IOobjectList& objects,
    const wordRes& selection,
    const bool filter
)
:
    classes_
    (
        filter
      ? objects.classes(selection)
      : objects.classes()
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::IOobjectNames::clear()
{
    classes_.clear();
}


bool Foam::IOobjectNames::empty() const
{
    forAllConstIters(classes_, iter)
    {
        if (!iter.object().empty())
        {
            return false;
        }
    }

    return true;
}


Foam::label Foam::IOobjectNames::size() const
{
    label total = 0;

    forAllConstIters(classes_, iter)
    {
        total += iter.object().size();
    }

    return total;
}


Foam::label Foam::IOobjectNames::count
(
    const UList<word>& whichTypes
) const
{
    return countImpl(classes_, whichTypes);
}


Foam::label Foam::IOobjectNames::count
(
    const wordHashSet& whichTypes
) const
{
    return countImpl(classes_, whichTypes);
}


void Foam::IOobjectNames::prune_0()
{
    forAllIters(classes_, iter)
    {
        iter.object().filterKeys
        (
            [](const word& k){ return k.endsWith("_0"); },
            true  // prune
        );
    }
}


void Foam::IOobjectNames::removeClasses(const UList<word>& whichTypes)
{
    classes_.erase(whichTypes);
}


void Foam::IOobjectNames::removeClasses(const wordHashSet& whichTypes)
{
    classes_.erase(whichTypes);
}


void Foam::IOobjectNames::retainClasses(const UList<word>& whichTypes)
{
    retainClassesImpl(whichTypes);
}


void Foam::IOobjectNames::retainClasses(const wordHashSet& whichTypes)
{
    retainClassesImpl(whichTypes);
}


bool Foam::IOobjectNames::checkNames() const
{
    if (Pstream::parRun())
    {
        // Count
        label len = 0;

        forAllConstIters(classes_, iter)
        {
            len += iter.object().size();
        }

        wordList masterNames(len);
        len = 0;

        forAllConstIters(classes_, iter)
        {
            const wordHashSet& objNames = iter.object();
            for (const word& objName : objNames)
            {
                masterNames[len] = objName;
                ++len;
            }
        }

        // Sort for consistent order
        Foam::sort(masterNames);

        const wordList localNames(masterNames);
        Pstream::scatter(masterNames);

        if (localNames != masterNames)
        {
            FatalErrorInFunction
                << "Objects not synchronised across processors." << nl
                << "Master has " << flatOutput(masterNames) << nl
                << "Processor " << Pstream::myProcNo()
                << " has " << flatOutput(localNames)
                << exit(FatalError);

            return false;
        }
    }

    return true;
}


void Foam::IOobjectNames::reduce()
{
    if (Pstream::parRun())
    {
        Pstream::mapCombineGather(classes_, HashSetOps::plusEqOp<word>());
        Pstream::mapCombineScatter(classes_);
    }
}


Foam::wordList Foam::IOobjectNames::sortedNames() const
{
    wordList list(this->size());

    label n = 0;

    forAllConstIters(classes_, iter)
    {
        for (const word& name : iter.object())
        {
            list[n] = name;
            ++n;
        }
    }

    inplaceUniqueSort(list);

    return list;
}


Foam::Ostream& Foam::IOobjectNames::info(Ostream& os) const
{
    for (const word& clsName : classes_.sortedToc())
    {
        os  << clsName << "  "
            << flatOutput(classes_[clsName].sortedToc()) << nl;
    }

    return os;
}


// ************************************************************************* //
