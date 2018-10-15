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

#include "PatchFunction1.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

template<class Type>
Foam::PatchFunction1<Type>::PatchFunction1
(
    const polyPatch& pp,
    const word& entryName
)
:
    refCount(),
    name_(entryName),
    patch_(pp)
{}


template<class Type>
Foam::PatchFunction1<Type>::PatchFunction1(const PatchFunction1<Type>& pf1)
:
    refCount(),
    name_(pf1.name_),
    patch_(pf1.patch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::word& Foam::PatchFunction1<Type>::name() const
{
    return name_;
}


template<class Type>
void Foam::PatchFunction1<Type>::convertTimeBase(const Time&)
{}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::PatchFunction1<Type>::value
(
    const scalar x
) const
{
    NotImplemented;

    return Field<Type>();
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::PatchFunction1<Type>::integrate
(
    const scalar x1,
    const scalar x2
) const
{
    NotImplemented;

    return Field<Type>();
}


template<class Type>
void Foam::PatchFunction1<Type>::autoMap(const FieldMapper& mapper)
{
    NotImplemented;
}


template<class Type>
void Foam::PatchFunction1<Type>::rmap
(
    const PatchFunction1<Type>& pf1,
    const labelList& addr
)
{
    NotImplemented;
}


template<class Type>
void Foam::PatchFunction1<Type>::writeData(Ostream& os) const
{
    os.writeKeyword(name_) << type();
}


// * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const PatchFunction1<Type>& pf1
)
{
    os.check(FUNCTION_NAME);

    os  << pf1.name_;
    pf1.writeData(os);

    return os;
}


// ************************************************************************* //
