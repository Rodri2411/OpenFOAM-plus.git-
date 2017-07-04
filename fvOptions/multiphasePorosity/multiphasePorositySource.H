/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
     \\/     M anipulation  | Copyright (C) 2017 IH-Cantabria
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
    Foam::fv::multiphasePorositySource

Group
    grpFvOptionsSources

Description

Usage

SourceFiles
    multiphasePorositySource.C

\*---------------------------------------------------------------------------*/

#ifndef multiphasePorositySource_H
#define multiphasePorositySource_H

#include "fvOption.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class porosityModel;

namespace fv
{


/*---------------------------------------------------------------------------*\
                   Class multiphasePorositySource Declaration
\*---------------------------------------------------------------------------*/

class multiphasePorositySource
:
    public option
{

protected:

    // Protected data

        // Coefficients per cell zone

            //- a coefficient
            scalarList aZone_;

            //- b coefficient
            scalarList bZone_;

            //- c coefficient
            scalarList cZone_;

            //- D50 median diameter
            scalarList D50Zone_;

            //- Porosity coefficient
            scalarList porosityZone_;

            //- Zone indices
            labelListList zoneIDs_;


        // Field properties

            //- Name of laminar kinematic viscosity; default = nu
            word nuName_;

            //- Name of porosity field
            word porosityName_;

            //- Porosity field
            //  Note:  typically created by the solver
            volScalarField& porosity_;


        //- Keulegan-Carpenter parameter; negative value to switch off
        //  Default = -1
        scalar kc_;


private:

    // Private Member Functions

        //- Set the porosity field
        void setPorosity();

        //- Disallow default bitwise copy construct
        multiphasePorositySource(const multiphasePorositySource&);

        //- Disallow default bitwise assignment
        void operator=(const multiphasePorositySource&);


public:

    //- Runtime type information
    TypeName("multiphasePorositySource");


    // Constructors

        //- Construct from components
        multiphasePorositySource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~multiphasePorositySource()
    {}


    // Member Functions

        // Add explicit and implicit contributions

            //- Add implicit contribution to phase momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldi
            );


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
