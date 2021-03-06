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

#include "systemCall.H"
#include "Time.H"
#include "dynamicCode.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(systemCall, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        systemCall,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::label Foam::functionObjects::systemCall::dispatch(const stringList& calls)
{
    if (calls.empty())
    {
        return 0;
    }

    label nCalls = 0;

    if (!masterOnly_ || Pstream::master())
    {
        for (const string& call : calls)
        {
            Foam::system(call); // Handles empty command as a successful no-op.
            ++nCalls;
        }
    }

    // MPI barrier
    if (masterOnly_)
    {
        Pstream::scatter(nCalls);
    }

    return nCalls;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::systemCall::systemCall
(
    const word& name,
    const Time&,
    const dictionary& dict
)
:
    functionObject(name),
    executeCalls_(),
    endCalls_(),
    writeCalls_(),
    masterOnly_(false)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::systemCall::read(const dictionary& dict)
{
    functionObject::read(dict);

    executeCalls_.clear();
    writeCalls_.clear();
    endCalls_.clear();

    dict.readIfPresent("executeCalls", executeCalls_);
    dict.readIfPresent("writeCalls", writeCalls_);
    dict.readIfPresent("endCalls", endCalls_);
    masterOnly_ = dict.lookupOrDefault("master", false);

    if (executeCalls_.empty() && endCalls_.empty() && writeCalls_.empty())
    {
        WarningInFunction
            << "no executeCalls, endCalls or writeCalls defined."
            << endl;
    }
    else if (!dynamicCode::allowSystemOperations)
    {
        FatalErrorInFunction
            << "Executing user-supplied system calls may not be enabled by "
            << "default due to potential " << nl
            << "security issues.  If you trust the case you can enable this "
            << "facility by " << nl
            << "adding to the InfoSwitches setting in the system controlDict:"
            << nl << nl
            << "    allowSystemOperations 1" << nl << nl
            << "The system controlDict is either" << nl << nl
            << "    ~/.OpenFOAM/$WM_PROJECT_VERSION/controlDict" << nl << nl
            << "or" << nl << nl
            << "    $WM_PROJECT_DIR/etc/controlDict" << nl << nl
            << exit(FatalError);
    }

    return true;
}


bool Foam::functionObjects::systemCall::execute()
{
    dispatch(executeCalls_);
    return true;
}


bool Foam::functionObjects::systemCall::write()
{
    dispatch(writeCalls_);
    return true;
}


bool Foam::functionObjects::systemCall::end()
{
    dispatch(endCalls_);
    return true;
}


// ************************************************************************* //
