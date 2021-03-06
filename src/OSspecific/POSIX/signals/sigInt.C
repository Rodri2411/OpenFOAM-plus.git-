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

#include "sigInt.H"
#include "error.H"
#include "JobInfo.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

struct sigaction Foam::sigInt::oldAction_;

bool Foam::sigInt::sigActive_ = false;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sigInt::sigHandler(int)
{
    // Reset old handling
    if (sigaction(SIGINT, &oldAction_, nullptr) < 0)
    {
        FatalErrorInFunction
            << "Cannot reset SIGINT trapping"
            << abort(FatalError);
    }

    jobInfo.signalEnd();        // Update jobInfo file
    raise(SIGINT);              // Throw signal (to old handler)
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sigInt::sigInt()
{
    set(false);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sigInt::~sigInt()
{
    unset(false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sigInt::set(bool)
{
    if (!sigActive_)
    {
        struct sigaction newAction;
        newAction.sa_handler = sigHandler;
        newAction.sa_flags = SA_NODEFER;
        sigemptyset(&newAction.sa_mask);
        if (sigaction(SIGINT, &newAction, &oldAction_) < 0)
        {
            FatalErrorInFunction
                << "Cannot call sigInt::set() more than once"
                << abort(FatalError);
        }
        sigActive_ = true;
    }
}


void Foam::sigInt::unset(bool)
{
    if (sigActive_)
    {
        if (sigaction(SIGINT, &oldAction_, nullptr) < 0)
        {
            FatalErrorInFunction
                << "Cannot set SIGINT trapping"
                << abort(FatalError);
        }
        sigActive_ = false;
    }
}


// ************************************************************************* //
