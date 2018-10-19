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

#include "sampledSurfaceProxy.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::sampledNull Foam::sampledSurfaceProxy::dummy_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledSurfaceProxy::sampledSurfaceProxy()
:
    surf_(std::ref<sampledSurface>(dummy_)),
    merged_(),
    state_(updateState::TOPO_CHANGE)
{}


Foam::sampledSurfaceProxy::sampledSurfaceProxy(sampledSurface& s)
:
    surf_(std::ref<sampledSurface>(s)),
    merged_(),
    state_(updateState::TOPO_CHANGE)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sampledSurfaceProxy::surface(sampledSurface& s)
{
    surf_ = std::ref<sampledSurface>(s);
    merged_.clear();
    state_ = updateState::TOPO_CHANGE;
}


const Foam::sampledSurface& Foam::sampledSurfaceProxy::surface() const
{
    return surf_.get();
}


Foam::sampledSurface& Foam::sampledSurfaceProxy::surface()
{
    return surf_.get();
}


void Foam::sampledSurfaceProxy::merge(scalar mergeDim)
{
    if (merged_.merge(surface(), mergeDim))
    {
        state_ = updateState::STAGED;
    }
}


// void Foam::sampledSurfaceProxy::clear()
// {
//     surface().clear();
//     merged_.clear();
//     state_ = updateState::TOPO_CHANGE;
// }


bool Foam::sampledSurfaceProxy::expire()
{
    if (surface().expire())
    {
        merged_.clear();
        state_ = updateState::TOPO_CHANGE;
        return true;
    }

    return false;
}


bool Foam::sampledSurfaceProxy::update(scalar mergeDim)
{
    if (surface().update())
    {
        merged_.merge(surface(), mergeDim);
        state_ = updateState::STAGED;
        return true;
    }

    return false;
}


bool Foam::sampledSurfaceProxy::staged() const
{
    return (state_ == updateState::STAGED);
}


bool Foam::sampledSurfaceProxy::handled()
{
    if (state_ == updateState::STAGED)
    {
        state_ = updateState::UNCHANGED;
    }

    return (state_ == updateState::UNCHANGED);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const sampledSurfaceProxy& proxy)
{
    return (os << proxy.surface());
}


// ************************************************************************* //
