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
    Foam::fv::mangrovesSource

Group
    grpFvOptionsSources

Description
    Applies sink terms to momentum, k and epsilon equations to approximate
    flow interaction with mangroves.

Usage
    Example usage
    \verbatim
    mangroves
    {
        type            mangrovesSource;

        // Optional field names
        U               U;
        k               k;
        epsilon         epsilon;

        regions
        {
            region1
            {
                cellZone        c0;
                a               0.01;
                N               560;
                Ckp        	    1;
                Cep         	3.5;
                Cm        	    1;
                Cd              1.52;
            }
        }
    }
    \endverbatim

SourceFiles
    mangrovesSource.C

\*---------------------------------------------------------------------------*/

#ifndef mangrovesSource_H
#define mangrovesSource_H

#include "fvOption.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{


/*---------------------------------------------------------------------------*\
                       Class mangrovesSource Declaration
\*---------------------------------------------------------------------------*/

class mangrovesSource
:
    public option
{

protected:

    // Protected data

        // Coefficients per cell zone

            //- Width of the vegetation element [m]
            scalarList aZone_;

            //- Number of plants per unit of area [1/m2]
            scalarList NZone_;

            //- Ckp
            scalarList CkpZone_;

            //- Cep
            scalarList CepZone_;

            //- Inertia coefficient
            scalarList CmZone_;

            //- Drag coefficient
            scalarList CdZone_;

            //- Zone indices
            labelListList zoneIDs_;


        // Field properties

            //- Name of U; default = U
            word UName_;

            //- Name of k; default = k
            word kName_;

            //- Name of epsilon; default = epsilon
            word epsilonName_;


        // Private Member Functions

            //- Return the k contribution
            tmp<volScalarField> kMangrove() const;

            //- Return the epsilon contribution
            tmp<volScalarField> epsilonMangrove() const;

            //- Return the drag force contribution
            tmp<volScalarField> UDragMangrove() const;

            //- Return the inertia force contribution
            tmp<volScalarField> UInertiaMangrove() const;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        mangrovesSource(const mangrovesSource&);

        //- Disallow default bitwise assignment
        void operator=(const mangrovesSource&);


public:

    //- Runtime type information
    TypeName("mangrovesSource");


    // Constructors

        //- Construct from components
        mangrovesSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~mangrovesSource()
    {}


    // Member Functions

        // Add explicit and implicit contributions

            //- Add implicit contribution to k/epsilon equations
            virtual void addSup
            (
                fvMatrix<scalar>& eqn,
                const label fieldi
            );

            //- Add implicit contribution to compressible k/epsilon equations
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const label fieldi
            );

            //- Add implicit contribution to phase k/epsilon equations
            virtual void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const label fieldi
            );

            //- Add implicit contribution to momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldi
            );

            //- Add implicit contribution to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldi
            );

            //- Add implicit contribution to phase momentum equation
            virtual void addSup
            (
                const volScalarField& alpha,
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
