/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                isoAdvector | Copyright (C) 2016-2017 DHI
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
    setAlphaField

Description
    Uses isoCutCell to create a volume fraction field from either a cylinder,
    a sphere or a plane.

    Original code supplied by Johan Roenby, DHI (2016)

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "isoCutFace.H"
#include "isoCutCell.H"
#include "Enum.H"
#include "mathematicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class shapeSelector
{
    public:

        enum class shapeType
        {
            PLANE,
            SPHERE,
            CYLINDER,
            SIN
        };

    static const Foam::Enum<shapeType> shapeTypeNames;
};


const Foam::Enum
<
    shapeSelector::shapeType
>
shapeSelector::shapeTypeNames
{
    { shapeSelector::shapeType::PLANE, "plane" },
    { shapeSelector::shapeType::SPHERE, "sphere" },
    { shapeSelector::shapeType::CYLINDER, "cylinder" },
    { shapeSelector::shapeType::SIN, "sin" },
};


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    Info<< "Reading setAlphaFieldDict\n" << endl;

    IOdictionary dict
    (
        IOobject
        (
            "setAlphaFieldDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const shapeSelector::shapeType surfType
    (
        shapeSelector::shapeTypeNames.read(dict.lookup("type"))
    );
    const vector centre(dict.lookup("centre"));
    const word fieldName(dict.lookup("field"));

    Info<< "Reading field " << fieldName << "\n" << endl;
    volScalarField alpha1
    (
        IOobject
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    scalar f0 = 0.0;
    scalarField f(mesh.points().size());

    Info<< "Processing type '" << shapeSelector::shapeTypeNames[surfType]
        << "'" << endl;

    switch (surfType)
    {
        case shapeSelector::shapeType::PLANE:
        {
            const vector direction(dict.lookup("direction"));

            f = -(mesh.points() - centre) & (direction/mag(direction));
            f0 = 0.0;
            break;
        }
        case shapeSelector::shapeType::SPHERE:
        {
            const scalar radius(readScalar(dict.lookup("radius")));

            f = -mag(mesh.points() - centre);
            f0 = -radius;
            break;
        }
        case shapeSelector::shapeType::CYLINDER:
        {
            const scalar radius(readScalar(dict.lookup("radius")));
            const vector direction(dict.lookup("direction"));

            f = -sqrt
            (
                sqr(mag(mesh.points() - centre))
              - sqr(mag((mesh.points() - centre) & direction))
            );
            f0 = -radius;
            break;
        }
        case shapeSelector::shapeType::SIN:
        {
            const scalar period(readScalar(dict.lookup("period")));
            const scalar amplitude(readScalar(dict.lookup("amplitude")));
            const vector up(dict.lookup("up"));
            const vector direction(dict.lookup("direction"));

            const scalarField xx
            (
                (mesh.points() - centre) & direction/mag(direction)
            );
            const scalarField zz((mesh.points() - centre) & up/mag(up));

            f = amplitude*Foam::sin(2*mathematical::pi*xx/period) - zz;
            f0 = 0;
            break;
        }
    }


    // Define function on mesh points and isovalue

    // Calculating alpha1 volScalarField from f = f0 isosurface
    isoCutCell icc(mesh, f);
    icc.volumeOfFluid(alpha1, f0);

    // Writing volScalarField alpha1
    ISstream::defaultPrecision(18);
    alpha1.write();

    Info<< nl << "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") - 1 = " << max(alpha1).value() - 1
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
