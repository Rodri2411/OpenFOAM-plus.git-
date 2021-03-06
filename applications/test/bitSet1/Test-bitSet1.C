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

Application
    Test-bitSet1

Description
    Test bitSet functionality

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "boolList.H"
#include "bitSet.H"
#include "HashSet.H"
#include "cpuTime.H"
#include <vector>
#include <unordered_set>

using namespace Foam;

inline Ostream& report
(
    const bitSet& bitset,
    bool showBits = false,
    bool debugOutput = false
)
{
    Info<< "size=" << bitset.size() << "/" << bitset.capacity()
        << " count=" << bitset.count()
        << " all:" << bitset.all()
        << " any:" << bitset.any()
        << " none:" << bitset.none() << nl;

    Info<< "values: " << flatOutput(bitset) << nl;
    if (showBits)
    {
        bitset.printBits(Info, debugOutput) << nl;
    }

    return Info;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
    bitSet set1(100);

    Info<<"bitSet(label): "; report(set1, true);

    bitSet set2(100, { -1, 10, 25, 45});
    Info<<"bitSet(label, labels): "; report(set2, true);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
