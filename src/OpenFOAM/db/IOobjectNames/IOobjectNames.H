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
    Foam::IOobjectNames

Description
    Simple management of IOobject names, usually arising from an IOobjectList.

SourceFiles
    IOobjectNames.C
    IOobjectNamesTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef IOobjectNames_H
#define IOobjectNames_H

#include "HashSet.H"
#include "wordRes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class objectRegistry;
class IOobjectList;

/*---------------------------------------------------------------------------*\
                        Class IOobjectNames Declaration
\*---------------------------------------------------------------------------*/

class IOobjectNames
{
    // Private Member Data

        //- IOobject names, categorized by classes
        HashTable<wordHashSet> classes_;


    // Private Member Functions

        //- Implementation for count()
        template<class Container>
        static label countImpl
        (
            const HashTable<wordHashSet>& classes,
            const Container& whichTypes
        );

        //- Implementation for retainClasses()
        template<class Container>
        void retainClassesImpl(const Container& whichTypes);


public:

    // Constructors

        //- Construct null
        IOobjectNames();

        //- Construct copy of objectRegistry::classes()
        explicit IOobjectNames(const objectRegistry& objects);

        //- Construct copy of IOobjectList::classes()
        explicit IOobjectNames(const IOobjectList& objects);

        //- Construct copy of objectRegistry::classes(),
        //- optionally filtered by object name
        IOobjectNames
        (
            const objectRegistry& objects,
            const wordRes& selection,
            const bool filter = true
        );

        //- Construct copy of IOobjectList::classes(),
        //- optionally filtered by object name
        IOobjectNames
        (
            const IOobjectList& objects,
            const wordRes& selection,
            const bool filter = true
        );


    //- Destructor
    ~IOobjectNames() = default;


    // Static Member Functions

        //- Count fields of known types matching the filter selection
        static label count
        (
            const objectRegistry& objects,
            const UList<word>& knownTypes,
            const wordRes& selection,
            const bool filter = true
        );

        //- Count fields of known types matching the filter selection
        static label count
        (
            const objectRegistry& objects,
            const wordHashSet& knownTypes,
            const wordRes& selection,
            const bool filter = true
        );

        //- Count fields of known types matching the filter selection
        static label count
        (
            const IOobjectList& objects,
            const UList<word>& knownTypes,
            const wordRes& selection,
            const bool filter = true
        );

        //- Count fields of known types matching the filter selection
        static label count
        (
            const IOobjectList& objects,
            const wordHashSet& knownTypes,
            const wordRes& selection,
            const bool filter = true
        );


    // Member Functions

        //- Clear all entries
        void clear();

        //- True if there are no object names
        bool empty() const;

        //- Total number of object names
        label size() const;


        //- Total number of object names of specified classes
        label count(const UList<word>& whichTypes) const;

        //- Total number of object names of specified classes
        label count(const wordHashSet& whichTypes) const;


        //- Filter fields based on their names, keeping or pruning entries
        //- that satisfy the unary predicate.
        //
        //  \return The number of items changed (removed)
        template<class UnaryPredicate>
        label filter(const UnaryPredicate& pred, bool pruning = false);

        //- Remove field names ending with "_0" (restart fields)
        void prune_0();

        //- Remove the specified classes
        void removeClasses(const UList<word>& whichTypes);

        //- Remove the specified classes
        void removeClasses(const wordHashSet& whichTypes);

        //- Retain the specified classes only
        void retainClasses(const UList<word>& whichTypes);

        //- Retain the specified classes only
        void retainClasses(const wordHashSet& whichTypes);


        //- Verify that names are synchronised across processors
        bool checkNames() const;

        //- Synchronised entries across processors
        //  Ensure all processes have identical information
        void reduce();


        //- The sorted names of IOobjects
        wordList sortedNames() const;

        //- The sorted names of IOobjects with the given ClassType
        template<class ClassType>
        wordList sortedNames() const;

        //- Information
        Ostream& info(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "IOobjectNamesTemplates.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
