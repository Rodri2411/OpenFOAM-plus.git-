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

\*---------------------------------------------------------------------------*/

#include "IOstream.H"
#include "axesRotation.H"
#include "coordinateSystem.H"
#include "coordinateSystems.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coordinateSystem, 0);
    defineRunTimeSelectionTable(coordinateSystem, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordinateSystem::coordinateSystem()
:
    name_(),
    note_(),
    origin_(Zero),
    R_(new axesRotation(sphericalTensor::I))
{}


Foam::coordinateSystem::coordinateSystem(const coordinateSystem& cs)
:
    name_(cs.name_),
    note_(cs.note_),
    origin_(cs.origin_),
    R_(cs.R_.clone())
{}


Foam::coordinateSystem::coordinateSystem(coordinateSystem&& cs)
:
    name_(std::move(cs.name_)),
    note_(std::move(cs.note_)),
    origin_(std::move(cs.origin_)),
    R_(std::move(cs.R_))
{}


Foam::coordinateSystem::coordinateSystem
(
    const word& name,
    const coordinateSystem& cs
)
:
    name_(name),
    note_(cs.note_),
    origin_(cs.origin_),
    R_(cs.R_.clone())
{}


Foam::coordinateSystem::coordinateSystem
(
    const word& name,
    const point& origin,
    const coordinateRotation& cr
)
:
    name_(name),
    note_(),
    origin_(origin),
    R_(cr.clone())
{}


Foam::coordinateSystem::coordinateSystem
(
    const word& name,
    const point& origin,
    const vector& axis,
    const vector& dirn
)
:
    name_(name),
    note_(),
    origin_(origin),
    R_(new axesRotation(axis, dirn))
{}


Foam::coordinateSystem::coordinateSystem
(
    const word& name,
    const dictionary& dict
)
:
    name_(name),
    note_(),
    origin_(Zero),
    R_()
{
    init(dict);
}


Foam::coordinateSystem::coordinateSystem(const dictionary& dict)
:
    name_(),
    note_(),
    origin_(Zero),
    R_()
{
    init(dict);
}


Foam::coordinateSystem::coordinateSystem
(
    const objectRegistry& obr,
    const dictionary& dict
)
:
    name_(),
    note_(),
    origin_(Zero),
    R_()
{
    const entry* entryPtr = dict.lookupEntryPtr(typeName_(), false, false);

    if (!entryPtr)
    {
        // No 'coordinateSystem' entry
        init(dict, obr);
    }
    else if (entryPtr->isDict())
    {
        // 'coordinateSystem' as dictionary entry - use it
        init(entryPtr->dict(), obr);
    }
    else
    {
        // 'coordinateSystem' as non-dictionary entry
        // - this is a lookup into global coordinateSystems

        keyType key(entryPtr->stream());

        const coordinateSystems& lst = coordinateSystems::New(obr);
        const label index = lst.findIndex(key);

        if (debug)
        {
            InfoInFunction
                << "Using global coordinate system: "
                << key << "=" << index << endl;
        }

        if (index < 0)
        {
            FatalErrorInFunction
                << "could not find coordinate system: " << key << nl
                << "available coordinate systems: " << lst.toc() << nl << nl
                << exit(FatalError);
        }

        // Copy from coordinateSystem, but assign the name as the typeName
        // to avoid strange things in writeDict()
        operator=(lst[index]);
        name_ = typeName_();
    }
}


Foam::coordinateSystem::coordinateSystem(Istream& is)
:
    name_(is),
    note_(),
    origin_(Zero),
    R_()
{
    dictionary dict(is);
    init(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dictionary Foam::coordinateSystem::dict(bool ignoreType) const
{
    dictionary dict;

    dict.add("name", name_);

    // Only write type for derived types
    if (!ignoreType && type() != typeName_())
    {
        dict.add("type", type());
    }

    // The note entry is optional
    if (note_.size())
    {
        dict.add("note", note_);
    }

    dict.add("origin", origin_);
    dict.add("e1", R_->e1());
    dict.add("e3", R_->e3());

    return dict;
}


Foam::vector Foam::coordinateSystem::localToGlobal
(
    const vector& local,
    bool translate
) const
{
    if (translate)
    {
        return (R_->transform(local)) + origin_;
    }

    return R_->transform(local);
}


Foam::tmp<Foam::vectorField> Foam::coordinateSystem::localToGlobal
(
    const vectorField& local,
    bool translate
) const
{
    if (translate)
    {
        return (R_->transform(local)) + origin_;
    }

    return R_->transform(local);
}


Foam::vector Foam::coordinateSystem::globalToLocal
(
    const vector& global,
    bool translate
) const
{
    if (translate)
    {
        return R_->invTransform(global - origin_);
    }

    return R_->invTransform(global);
}


Foam::tmp<Foam::vectorField> Foam::coordinateSystem::globalToLocal
(
    const vectorField& global,
    bool translate
) const
{
    if (translate)
    {
        return R_->invTransform(global - origin_);
    }

    return R_->invTransform(global);
}


void Foam::coordinateSystem::clear()
{
    note_.clear();
    origin_ = Zero;
    R_->clear();
}


void Foam::coordinateSystem::transfer(coordinateSystem& cs)
{
    name_ = std::move(cs.name_);
    note_ = std::move(cs.note_);
    origin_ = std::move(cs.origin_);
    R_ = std::move(cs.R_);
}


void Foam::coordinateSystem::write(Ostream& os) const
{
    os  << type() << " origin: " << origin() << nl;
    R_->write(os);
}


void Foam::coordinateSystem::writeDict(Ostream& os, bool subDict) const
{
    if (subDict)
    {
        os.beginBlock(name_);
    }

    os.writeEntry("type", type());

    if (note_.size())
    {
        // The 'note' is optional
        os.writeEntry("note", note_);
    }

    os.writeEntry("origin", origin_);
    R_->write(os);

    if (subDict)
    {
        os.endBlock();
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::coordinateSystem::operator=(const coordinateSystem& cs)
{
    name_ = cs.name_;
    note_ = cs.note_;
    origin_ = cs.origin_;

    // Some extra safety
    if (cs.R_.valid())
    {
        R_ = cs.R_.clone();
    }
    else
    {
        R_.reset(new axesRotation(sphericalTensor::I));
    }
}

void Foam::coordinateSystem::operator=(coordinateSystem&& cs)
{
    transfer(cs);
}


void Foam::coordinateSystem::init(const dictionary& dict)
{
    dict.lookup("origin") >> origin_;
    note_.clear();
    dict.readIfPresent("note", note_);
    R_ = coordinateRotation::New(dict.subDict("coordinateRotation"));
}


void Foam::coordinateSystem::init
(
    const dictionary& dict,
    const objectRegistry& obr
)
{
    dict.lookup("origin") >> origin_;

    // The 'note' entry is optional
    note_.clear();
    dict.readIfPresent("note", note_);
    R_ = coordinateRotation::New(dict.subDict("coordinateRotation"), obr);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

bool Foam::operator!=(const coordinateSystem& a, const coordinateSystem& b)
{
    return
    (
        a.origin() != b.origin()
     || a.type() != b.type()
     || a.R().R() != b.R().R()
    );
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const coordinateSystem& cs)
{
    cs.write(os);
    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
