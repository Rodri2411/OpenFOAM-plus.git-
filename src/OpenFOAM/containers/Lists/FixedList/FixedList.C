/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

#include "FixedList.H"
#include "ListLoopM.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, unsigned Size>
Foam::label Foam::FixedList<T, Size>::find
(
    const T& val,
    const label start
) const
{
    if (start >= 0)
    {
        List_CONST_ACCESS(T, *this, list);

        for (label i = start; i < label(Size); ++i)
        {
            if (list[i] == val)
            {
                return i;
            }
        }
    }

    return -1;
}


template<class T, unsigned Size>
void Foam::FixedList<T, Size>::moveFirst(const label i)
{
    checkIndex(i);

    for (label lower = 0; lower < i; ++lower)
    {
        Foam::Swap(v_[lower], v_[i]);
    }
}


template<class T, unsigned Size>
void Foam::FixedList<T, Size>::moveLast(const label i)
{
    checkIndex(i);

    for (label upper = label(Size - 1); upper > i; --upper)
    {
        Foam::Swap(v_[i], v_[upper]);
    }
}


template<class T, unsigned Size>
void Foam::FixedList<T, Size>::swapFirst(const label i)
{
    checkIndex(i);

    if (i > 0)
    {
        Foam::Swap(v_[0], v_[i]);
    }
}


template<class T, unsigned Size>
void Foam::FixedList<T, Size>::swapLast(const label i)
{
    checkIndex(i);

    const label upper = label(Size - 1);

    if (i < upper)
    {
        Foam::Swap(v_[i], v_[upper]);
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, unsigned Size>
bool Foam::FixedList<T, Size>::operator==(const FixedList<T, Size>& list) const
{
    List_CONST_ACCESS(T, *this, lhs);
    List_CONST_ACCESS(T, (list), rhs);

    // List sizes are identical by definition (template parameter)
    for (unsigned i = 0; i < Size; ++i)
    {
        if (!(lhs[i] == rhs[i]))
        {
            return false;
        }
    }

    // Contents appear to be identical.
    return true;
}


template<class T, unsigned Size>
bool Foam::FixedList<T, Size>::operator<(const FixedList<T, Size>& list) const
{
    List_CONST_ACCESS(T, *this, lhs);
    List_CONST_ACCESS(T, (list), rhs);

    // List sizes are identical by definition (template parameter)
    for (unsigned i=0; i<Size; ++i)
    {
        if (lhs[i] < rhs[i])
        {
            return true;
        }
        else if (rhs[i] < lhs[i])
        {
            return false;
        }
    }

    // Contents appear to be identical.
    return false;
}


template<class T, unsigned Size>
bool Foam::FixedList<T, Size>::operator!=(const FixedList<T, Size>& list) const
{
    return !operator==(list);
}


template<class T, unsigned Size>
bool Foam::FixedList<T, Size>::operator>(const FixedList<T, Size>& list) const
{
    return list.operator<(*this);
}


template<class T, unsigned Size>
bool Foam::FixedList<T, Size>::operator<=(const FixedList<T, Size>& list) const
{
    return !list.operator<(*this);
}


template<class T, unsigned Size>
bool Foam::FixedList<T, Size>::operator>=(const FixedList<T, Size>& list) const
{
    return !operator<(list);
}


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "FixedListIO.C"

// ************************************************************************* //
